classdef VRClass
   methods(Static)
        
       %This function finds the Eucledian distances
       function distance = findEuclideanDistance(one, two)
            [rowOne, colOne] = size(one);
            [rowTwo, colTwo] = size(two);
            
            if (rowOne == rowTwo)
                distance = zeros(colOne, colTwo);
                for ii=1:colOne
                    for jj=1:colTwo
                        %calculate the distances between each element of
                        %the arrays
                        distance(ii,jj)=sum((one(:,ii)-two(:,jj)).^2).^0.5; 
                    end
                end
            end
       end
       
       %This function computes the MFCC coefficients
       function result = mfcc(inputData, fs)
            distanceBtwBeginning = 100;    %distance between the beginning of each frame
            samplesPerFrame = 256;         %number of samples per frame
            numberOfFilters = 20;          %number of filters in the filter bank
            
            
            %split the audio file into frames
            [rows, cols] = size(inputData);
            numberOfFrames = floor((rows - samplesPerFrame) / distanceBtwBeginning) + 1;
            
            %getting samples of the input data
            for i = 1:samplesPerFrame
                for j = 1:numberOfFrames
                    samples(i, j) = inputData(((j - 1) * distanceBtwBeginning) + i); 
                end
            end 
            
            %applying the hamming window function
            HamSamples = diag(hamming(samplesPerFrame)) * samples; 
            
            %getting the spectrum 
            for i = 1:numberOfFrames
                spectrum(:, i) = fft(HamSamples(:, i)); %compute FFTs
            end
            
            
            % Find the mel filter bank values
            
            % We took the melfb function from Minh N. Do because it was 
            % too difficult to figure out and used it as it explained
            % in the comments to the function
            melAmplitudes = melfb(numberOfFilters, samplesPerFrame, fs);
            n2 = 1 + floor(samplesPerFrame / 2);
            melSpectrum = melAmplitudes * abs(spectrum(1:n2, :)).^2;
            
            % Get Discrete Cosine Transform to get MFCCs in order to
            % convert the spectrum back to the time domain
            result = dct(log(melSpectrum)); 
       end
       
       
       
       %This function assembles the codebook by vector quantization
       function codebook = makeCodebook(vectors,numCentroids)
           
            splitParameter = 0.03;          %splitting parameter
            codebook = mean(vectors, 2);    %initializing the codebook
            updatedDistortion = 10000;      %the initial value for the distortion  
                                            %in a cluster
            
            for i = 1:log2(numCentroids)
                %split/double the codebook based on the LBG rule
                codebook = [codebook*(1+splitParameter), codebook*(1-splitParameter)]; 
                
                distortion = 0;
                
                % while distortion is not too big then update keep looking for
                % centroid making the cluster
                while (updatedDistortion - distortion)/distortion > splitParameter
                    %find the distance between each codeword and quatization
                    %point
                    distance = VRClass.findEuclideanDistance(vectors, codebook);
                    
                    %find closest centroid to each training vector
                    [other, close] = min(distance, [], 2); 
                                                    
                    
                    
                    for j = 1:2^i
                        %updating the centroids
                        
                        codebook(:, j) = mean(vectors(:, find(close == j)), 2); 
                        
                        cluster = VRClass.findEuclideanDistance(vectors(:, find(j == close)), codebook(:, j)); 
                        [clusterRows, clusterCols] = size(cluster);
                        
                        %getting the distortion
                        for k = 1:clusterRows
                            distortion = distortion + cluster(k);
                        end
                        
                    end
                    %updating the distortion
                    updatedDistortion = distortion;
                    
                end
            end
            
       end
       
       
       %This function computes the distances and returns the closest 
       %audio file file from the library if there is a match within the threshold 
       function computeAndPrint(inputFileName, arrOfLibFiles, arrOfLibLabels, threshold, numCentroids)
            
            %getting the data from the input file
            [inputFileData, Fs] = audioread(inputFileName);
            %getting the mfcc from the input file
            inputFileMfcc = VRClass.mfcc(inputFileData,Fs);
            %getting the vector quatization coefficients from the input
            %file
            inputFileVQ = VRClass.makeCodebook(inputFileMfcc,numCentroids);
            
            %length of the library array
            [rows, length] = size(arrOfLibFiles);
                
            minDist = Inf;
            distanceIndex = 0;
            for i = 1:length
                %reading the files in the library and computing the MFCCs
                [audioData{i}, Fs] = audioread(arrOfLibFiles{i}); 
                mfccData{i} = VRClass.mfcc(audioData{i},Fs);
                
                %making the codebook data for each library vector
                vqData{i} = VRClass.makeCodebook(mfccData{i},numCentroids);
                
                %finding the distance between the input file and the
                %library files
                distanceMatrix{i} = VRClass.findEuclideanDistance(vqData{i}, inputFileVQ);
                
                %computing the distance matrix
                distances{i} = sum(min(distanceMatrix{i},[],2)) / size(distanceMatrix{i},1);
                
                
                %find minimum distance
                if distances{i} < minDist
                    minDist = distances{i};
                    distanceIndex = i;
                end
            end
            
            %display the result
            if minDist < threshold
                display('The closest sound from the library is:');
                display(arrOfLibFiles{distanceIndex});
                display('The closest label:');
                display(arrOfLibLabels{distanceIndex});
            else
                display('There is no match within the threshold');
            end

       end
       
   end    
end    