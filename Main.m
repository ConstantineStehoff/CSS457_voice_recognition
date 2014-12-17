clear all;
queryFileName = 'hello.wav';
fileNames = {'guitar.wav', 'hello.wav', 'meow.wav', 'hiss.wav', 'on.wav'};
labels = {'guitar', 'speech', 'cat meow', 'hiss', 'on'};
threshold = 10;
numberOfCentroids = 18;

VRClass.computeAndPrint(queryFileName, fileNames, labels, threshold, numberOfCentroids);


