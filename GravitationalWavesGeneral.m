%Band pass filter programme.

%Setting the data directory for imported data.
SimulationName='BU1GRNonRadialStergl2m2HLLE';
DataDirectory = strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data');
WebDirectory = strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'WebPageFile');

%Setting The Filter Characteristics
Peak = '1'; %1,2,3,..
FreqLow = 1850; FreqHigh = 2150;
FreqLowSt = num2str(FreqLow);FreqHighSt = num2str(FreqHigh);
FourierRangeInitial=300;FourierRangeFinal=10000;

%Important Constants In CGS Units.
GravConstCGS = 6.6738*(10^(-8)); %cm^3.g^-1.s^-2
LightSpeedCGS = 29979245800; %cm.s^-1
SolarMassCGS = 1.9891*(10^33); %g

%Solar Masses To Seconds And Vice Versa Transformation Coefficients.
SolMasToSec = 4.92686*(10^(-6));
SecToSolMas = 1/SolMasToSec;

%Solar Masses To Meters And Vice Versa Transformation Coefficients.
SolMasToMeters = 1477.04;
MetersToSolMas = 1/SolMasToMeters;

for L=2:2
    for M=-L:L
        for R=130:10:130
            for z=2:3
                
               %Specifying data to be imported and the frequency window in Hertz.
               l=num2str(L);m=num2str(M);r=num2str(R); %Psi4 Indices.
    
               %Importing data from data directory.
               if z == 1
                  Type='General';
                  Psi4TimeSeriesReal = dlmread(strcat(DataDirectory,'/','YlmFiles/Psi4lmrRealPartListSecASD',l,m,'r',r,'.00.asc'));
                  Psi4TimeSeriesImaginary = dlmread(strcat(DataDirectory,'/','YlmFiles/Psi4lmrImaginaryPartListSecASD',l,m,'r',r,'.00.asc'));
                  IntDataListTimeSec=Psi4TimeSeriesReal(:,1);
                  IntDataListReal=Psi4TimeSeriesReal(:,2);
                  IntDataListImaginary=Psi4TimeSeriesImaginary(:,2);
                  IntDataListComplex=IntDataListReal+IntDataListImaginary*i;
                  Psi4TimeSeries=[IntDataListTimeSec IntDataListComplex];
               elseif z==2
                  Type='Real';
                  Psi4TimeSeries = dlmread(strcat(DataDirectory,'/','YlmFiles/Psi4lmrRealPartListSecASD',l,m,'r',r,'.00.asc'));
               elseif z==3
                  Type='Imaginary';
                  Psi4TimeSeries = dlmread(strcat(DataDirectory,'/','YlmFiles/Psi4lmrImaginaryPartListSecASD',l,m,'r',r,'.00.asc'));
               end

               %Total data points.
               TotalDataPoints = numel(Psi4TimeSeries(:,1));

               %Data points list.
               DataList = Psi4TimeSeries(:,2);

               %Time Lists is seconds and solar masses.
               TimeListSec = Psi4TimeSeries(:,1);
               TimeListSM = TimeListSec*SecToSolMas;

               %Time step is seconds and solar masses.
               TimeStepSec = TimeListSec(2);
               TimeStepSM = TimeStepSec*SecToSolMas;

               %Fourier Constants.
               FourierSampRate = 1/TimeStepSec;
               FourierHertzConv = FourierSampRate/TotalDataPoints;
               FourierMaxResFreqDU = TotalDataPoints/2;
               FourierMaxResFreqHz = TotalDataPoints/(2*TimeStepSec);
               FourierFreqResolution = round(1/(TimeStepSec*TotalDataPoints));

               %Fourier lists.
               AbsFourierData = 2*abs(fft(DataList))/TotalDataPoints;
               FouRangeInit = round(FourierRangeInitial/FourierHertzConv); FouRangeFin = round(FourierRangeFinal/FourierHertzConv);
               MaxAmpFou = 0;MaxIndFou = 0;
               for k=FouRangeInit:FouRangeFin
                   if AbsFourierData(k) > MaxAmpFou
                      MaxAmpFou = AbsFourierData(k);
                      MaxIndFou = k;
                   end
               end
               FrequencyListDU = (0:1:TotalDataPoints-1);
               FrequencyListHertz = FrequencyListDU*FourierHertzConv;

               %Band Pass Filter.
               FilterOrder = 2;
               [ButtB,ButtA] = butter(FilterOrder, [FreqLow FreqHigh]/(FourierSampRate/2),'bandpass');
               FilteredData = filter(ButtB,ButtA,DataList); 
              
               %Fourier lists (Filtered data).
               AbsFourierFilteredData =  2*abs(fft(FilteredData))/TotalDataPoints;
               FouRangeInitFilt=round(FourierRangeInitial/FourierHertzConv) ;FouRangeFinFilt=round(FourierRangeFinal/FourierHertzConv);
               MaxAmpFouFilt=0;MaxIndFouFilt=0;
               for k=FouRangeInitFilt:FouRangeFinFilt
                   if AbsFourierFilteredData(k) > MaxAmpFouFilt
                      MaxAmpFouFilt = AbsFourierFilteredData(k);
                      MaxIndFouFilt = k;
                   end
               end
               dlmwrite(strcat(DataDirectory,'/YlmFilteredData/YlmFilteredData',Type,'Yl',l,'m',m,'r',r,'Peak',Peak,'.asc'),[TimeListSec FilteredData],'delimiter','\t');

               %Graphs data vs time in seconds and corresponding FFTs in Hertz.
               GraphFouRange = FourierRangeFinal;
               GraphFouRangeFilt = GraphFouRange;

               figure;

               subplot(2,2,1)
               plot(TimeListSec,DataList)
               title('Unfiltered data');
               xlabel('t(s)');
               ylabel('\psi4');

               subplot(2,2,2)
               plot(FrequencyListHertz,AbsFourierData);
               title('FFT of unfiltered data');
               xlabel('f(Hz)');
               ylabel('fft values');
               axis([0 GraphFouRange 0 MaxAmpFou+0.1*MaxAmpFou])

               subplot(2,2,3)
               plot(TimeListSec,FilteredData);
               title('Filtered data');
               xlabel('t(s)');
               ylabel('\psi4');

               subplot(2,2,4)
               plot(FrequencyListHertz,AbsFourierFilteredData);
               title('FFT of filtered data');
               xlabel('f(Hz)');
               ylabel('fft values');
               axis([0 GraphFouRangeFilt 0 MaxAmpFouFilt+0.1*MaxAmpFouFilt])

               suptitle({strcat(Type,' part of \psi4');strcat('Distance of detector from source r=',r,' solar masses. Indices l=',l,' m=',m);strcat('Frequency window:',FreqLowSt,'-',FreqHighSt,'Hz')});
               print('-dpng',strcat(WebDirectory,'/YlmFiltered/FiltDataGraph',Type,'Yl',l,'m',m,'r',r,'Peak',Peak,'.png'));
              
            end
        end
    end
end
%freqz(ButtB, ButtA)