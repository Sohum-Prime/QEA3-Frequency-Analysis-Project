% Chordination
%
% A fast algorithm which reads audio files and determines the chords
% being played. This algorithm utilizes the Pitch Class Profile (PCP)
% algorithm developed by Takuya Fujishima ("Realtime Chord Recognition
% of Musical Sound: A System Using Cacommon Lisp Music", 1999).
%
%
% Syntax:
% detectedChords = Chordination(audioFile);
% detectedChords = Chordination(audioFile, plotDisplay);
% detectedChords = Chordination(audioFile, plotDisplay, liveDisplay);
% detectedChords = Chordination(audioFile, plotDisplay, liveDisplay, ...
%                               secondsPerWindow);
% detectedChords = Chordination(audioFile, plotDisplay, liveDisplay, ...
%                               secondsPerWindow, numOverlaps);
% detectedChords = Chordination(audioFile, plotDisplay, liveDisplay, ...
%                               secondsPerWindow, numOverlaps, ...
%                               dictionaryType);
% 
% Parameters:
% audioFile: the filename of a .mp3 or .wav file containing music
% secondsPerWindow: length of the windows used (defaults to .25 seconds)
% numOverlaps: number of overlapping windows (defaults to 10)
% liveDisplay: Boolean variable: if set to 1, plays the song and outputs
%              the detected chords to the command window (defaults to 0)
% plotDisplay: Boolean variable: if set to 1, outputs a plot displaying the
%              notes detected in each window (defaults to 0)
% dictionaryType: A string variable which is either 'simple' or 'full'; if
%                 'simple', uses a small dictionary of 64 chords; if
%                 'full', uses a dictionary of 660 chords.
%
% Output:
% detectedChords: a string vector containing the names of the chord
%                 detected in each window
% Ashok Aggarwal, Miheer Patankar, Rohan Paul, Brian Hughes, Phil Sisk
% EECS 351
% December 11, 2017

function detectedChords = Chordination(audioFile, plotDisplay, ...
                                       liveDisplay, secondsPerWindow, ...
                                       numOverlaps, dictionaryType)
    
    % Default values
    defaultSecondsPerWindow = .25;
    defaultNumOverlaps = 10;

    if nargin > 6 || nargin < 1
        error("Invalid number of arguments.");
    end
    switch nargin
        case 1
            secondsPerWindow = defaultSecondsPerWindow;
            numOverlaps = defaultNumOverlaps;
            liveDisplay = 0;
            plotDisplay = 0;
            dictionaryType = 'simple';
        case 2
            secondsPerWindow = defaultSecondsPerWindow;
            numOverlaps = defaultNumOverlaps;
            liveDisplay = 0;
            dictionaryType = 'simple';
        case 3
            secondsPerWindow=.25;
            numOverlaps = 10;
            dictionaryType = 'simple';
        case 4
            numOverlaps = defaultNumOverlaps;
            dictionaryType = 'simple';
        case 5
            dictionaryType = 'simple';
    end

    if strcmp(dictionaryType, 'simple')
        s = load("chordMapSimple.mat");
        chordMap = s.chordMap;
    elseif strcmp(dictionaryType, 'full')
        s = load("chordMapFull.mat");
        chordMap = s.chordMap;
    else
        error("Invalid dictionary type.");
    end

    threshold = .01;                   % sensitivity to notes
    HPSweight = .5;                    % weight of HPS being used

    % Smoothing filter parameters
    smoothingYes = 1;                  % Whether smoothing is being used         
    smoothingLength = 3*numOverlaps;   % Length of smoothing vector

    % Bandpass filter parameters
    f_low = 40;    %Hz
    f_high = 1500; %Hz;

    [song,Fs] = audioread(audioFile);
    % In case of stereo track, change to mono
    if length(song(1, :)) == 2
        song = song(:, 1) + song(:, 2);
    end

    % Create hamming window
    windowlength=floor(secondsPerWindow*Fs);
    hammingwindow=hamming(windowlength); 

    filt = BandpassFilter(f_low, f_high, windowlength, Fs);

    % Pitch Class Profile Prep stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    l=1:windowlength/2;
    frequencies = Fs/windowlength*l;
    FrequencyBins = mod(round(12*log2(frequencies/110)),12);
    notes = 440*2.^((0:27)/12);

    for i = 1:length(FrequencyBins)
       notesDist = frequencies(i) - notes + 10.1;
       posDist = notesDist( find(notesDist >= 0) );
       closeNote = min(posDist);

       if closeNote > 20
           FrequencyBins(i) = -1;
       end

    end
    FrequencyBins(1)=-1; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    pcpmatrix=zeros(12,floor(length(song)/windowlength*numOverlaps));
    mostLikelyChordKeys = zeros(1, length(pcpmatrix));
    confidenceLevels  = mostLikelyChordKeys;
    detectedChords = string(zeros(1, length(pcpmatrix)));

    % Perform chord detection on each windowed audio sample
    n = 1;
    for j=0:(windowlength/numOverlaps):length(song)-windowlength
        %disp(['Window ', num2str(j/windowlength*N)]);
        % Window sample

        y=hammingwindow.*song((floor(j)+1):(floor(j)+windowlength));    
        ffty=fft(y)/max(abs(fft(y)));

        % Bandpass filter
        ffty = filt.*ffty;

        % Find Harmonic Product Spectrum (HPS)
        ffty = HPSweight*HarmonicProductSpectrumFullLength(ffty, 2)...
            + (1 - HPSweight)*ffty;

        % Find Pitch Class Profile (PCP)
        pcp = SelectivePitchClassProfile(ffty, FrequencyBins);

        pcpmatrix(:,floor(j/windowlength*numOverlaps+1))=pcp/max(pcp);
        n = n + 1;
    end

    % Smooth the detected notes
    if smoothingYes
        smoothingFilter = zeros (1, smoothingLength);
        smoothingFilter((ceil(smoothingLength/2)+ 1):smoothingLength) = ...
            ones(1, floor(smoothingLength/2));
        for n = 1:12
            pcpmatrix(n, :) = conv(pcpmatrix(n, :),...
                smoothingFilter, 'same');
        end
    end

    % Determine most likely chord based on PCP
    for n = 1:length(pcpmatrix)
       [detectedChords(n), pcpmatrix(:, n), confidenceLevels(n)] =...
           ChordFinder(pcpmatrix(:, n), chordMap, threshold);
    end

     
    % Play the song and display the corresponding detected chords
    if liveDisplay
        sound(song, Fs);
        for i = 1:length(detectedChords)
            t = ["Detected Chord: " detectedChords(i)];
            disp(t);
            pause(secondsPerWindow/numOverlaps);
        end
    end

    % Plot the detected notes versus time
    if plotDisplay
        figure; 
        pcpmatrix(pcpmatrix < threshold) = 0;
        xvector=repmat(1:12,length(pcpmatrix),1);
        timevector=repmat((1:length(pcpmatrix))',1,12)*secondsPerWindow/numOverlaps;
        a = ["A", "A#", "B","C","C#","D","D#","E","F","F#","G","G#"];
        pcpmatrix(pcpmatrix == 0) = NaN;
        plot3(timevector, xvector, pcpmatrix, 'LineWidth',...
            1.5, 'MarkerEdgeColor', 'r');
        set(gca,'Ytick',1:12,'yticklabel',a')
        xlabel('Time (s)');
        xlim([0, length(pcpmatrix)*secondsPerWindow/numOverlaps]);
        ylabel('Detected Notes');
        grid on
        view(-20, 86)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Function definitions                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% BandpassFilter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filter = BandpassFilter(f_low, f_high, len, Fs)

spacing = len/Fs;
filter = zeros(len, 1);
filter(ceil((len - floor(f_high*spacing))):...
    (floor((len - floor(f_low*spacing) -1))))...
    = ones(1, floor((f_high - f_low)*spacing));
filter(ceil(floor(f_low*spacing) + 1):(floor(f_high*spacing)))...
    = ones(1, floor((f_high - f_low)*spacing));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% HarmonicProductSpectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hps= HarmonicProductSpectrumFullLength(ffty, harmonics)

    %harmonics is the number of harmonics to be removed by HPS
    hps=zeros(length(ffty),1); 

     for i=1:floor(length(ffty)/2)

         hps(i)=abs(ffty(i));
         numharm=floor(log2(length(ffty)/2/i));

         if numharm>harmonics
             numharm=harmonics;
         end

         for j=1:numharm
             hps(i)=hps(i)*abs(ffty(i*2^j));
         end
     end
     
    hps = hps/max(hps);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% PitchClassProfile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pcp = SelectivePitchClassProfile(hps, FrequencyBins)
    
    pcp=zeros(12,1); 
      
    for i= 0:11
        pcp(i+1)=sum((norm(hps(FrequencyBins==i))).^2); 
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% ChordFinder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chordName, notes, confidenceLevel] = ChordFinder(notes,...
                                               chordMap, threshold)
    notes = notes';
    if sum(notes) > 0    
        notes = notes/max(notes);
    end
    
    notes(notes<=threshold)=0; %filter out by threshold val

    % Reward notes for being above threshold, but dont set them all to 1
    notes = (notes + .5*ceil(notes))/1.5;

    k = cell2mat(keys(chordMap));
    
    k1 = de2bi(k, 12, 'left-msb');
     
    diffs = zeros(length(k), 1);
    for i = 1:length(k)
        diffs(i) = sum(abs(k1(i, :) - notes));
    end
        
    [diffsSorted, diffInds] = sort(diffs, 'ascend');
    
    mostLikelyKeys = zeros(1, 3);
    
    mostLikelyKeys(:, 1:3) = k(diffInds(1:3));
    
    confidenceLevel = (12 - diffsSorted(1))/12*100;
    
    chordName = chordMap(mostLikelyKeys(1)); 
    
    notes = notes';
end
