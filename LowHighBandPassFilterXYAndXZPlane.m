%Band pass filter for XY and XZ plane.

%Setting the data directories for imported data.
SimulationName = 'BU0CowlNonRadiall2m2StergHLLEFull';
DataStorageDirectory = strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/');
StellarCharacteristicsDataDirectory = strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/StellarCharacteristics/');
SimulationSpecifics2dDataDirectory = strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_SimulationSpecifics/');

%Important constants in CGS units.
GravConstCGS = 6.6738*(10^(-8)); %cm^3.g^-1.s^-2
LightSpeedCGS = 29979245800; %cm.s^-1
SolarMassCGS = 1.9891*(10^33); %g

%Solar masses to seconds and vice versa transformation coefficients.
SolMasToSec = 4.92686*(10^(-6));
SecToSolMas = 1/SolMasToSec;

%Solar masses to meters and vice versa transformation coefficients.
SolMasToMeters = 1477.04;
MetersToSolMas = 1/SolMasToMeters;

%Setting the stellar characteristics.
StellarCharacteristicsList = dlmread(strcat(StellarCharacteristicsDataDirectory,'StellarCharacteristicsList.asc'));

MassNS = StellarCharacteristicsList(1,1);
EqRadiusNS =StellarCharacteristicsList(2,1);
StellarOblatnessNS = StellarCharacteristicsList(3,1);
PolRadiusNS = StellarCharacteristicsList(4,1);
OmegaNS = StellarCharacteristicsList(5,1);
RotFreqNS = StellarCharacteristicsList(6,1);
AngVelNS = StellarCharacteristicsList(7,1);
Kapa = StellarCharacteristicsList(8,1);
Gama = StellarCharacteristicsList(9,1);
RhoCentInit = StellarCharacteristicsList(10,1);
RhoAtmInit = StellarCharacteristicsList(11,1);

%Setting the simulation specifics grid characteristics.
SimulationSpecifics2dList = dlmread(strcat(SimulationSpecifics2dDataDirectory,'SimulationSpecifics2dList.asc'));

ItSampleStep2d = SimulationSpecifics2dList(1,1);
TotalDataPointsFull2d = SimulationSpecifics2dList(2,1);
TotalDataPoints2d = SimulationSpecifics2dList(3,1);
FinalItSampleStep2d = SimulationSpecifics2dList(4,1);
ItTimeStepSolMas2d = SimulationSpecifics2dList(5,1);
ItTimeStepSec2d = SimulationSpecifics2dList(6,1);
SampTimeStepSolMas2d = SimulationSpecifics2dList(7,1);
SampTimeStepSec2d = SimulationSpecifics2dList(8,1);
TotSimTimeInSolMas2d = SimulationSpecifics2dList(9,1);
TotSimTimeInSec2d = SimulationSpecifics2dList(10,1);

%Setting the grid origin.
OriginCoordinetsList = dlmread(strcat(SimulationSpecifics2dDataDirectory,'OriginCoordinetsList.asc'));

CentX = OriginCoordinetsList(1,1);
CentY = OriginCoordinetsList(2,1);
CentZ = OriginCoordinetsList(3,1);

%Setting the FFT characteristics.
FFTSpecifics2dList = dlmread(strcat(SimulationSpecifics2dDataDirectory,'FFTSpecifics2dList.asc'));

TimeListSec = dlmread(strcat(SimulationSpecifics2dDataDirectory,'TimeListSec2d.asc'));
TimeListSM = TimeListSec*SecToSolMas;

TotalDataPoints = numel(TimeListSec);
TimeStepSec = FFTSpecifics2dList(3,1);
TimeStepSM = FFTSpecifics2dList(2,1);        

FourierSampRate = 1/TimeStepSec;
FourierHertzConv = FourierSampRate/TotalDataPoints;
FourierMaxResFreqDU = TotalDataPoints/2;
FourierMaxResFreqHz = TotalDataPoints/(2*TimeStepSec);
FourierFreqResolution = round(1/(TimeStepSec*TotalDataPoints));

%Setting the filter characteristics and implementing it.
FilterType = 3; % 1 for Highpass, 2 for Lowpass, 3 for Bandpass
Peak = '2'; %1,2,3,..
FreqLow = 3900; FreqHigh = 4300;
FreqLowSt = num2str(FreqLow);FreqHighSt = num2str(FreqHigh);
FourierRangeInitial=300;FourierRangeFinal=20000;
if FilterType == 1 
  FilterOrder = 8;
  [ButtB,ButtA] = butter(FilterOrder, FreqHigh/(FourierSampRate/2),'high');
elseif FilterType == 2
  FilterOrder = 8;
  [ButtB,ButtA] = butter(FilterOrder, FreqHigh/(FourierSampRate/2));
elseif FilterType == 3
  FilterOrder = 2;
  [ButtB,ButtA] = butter(FilterOrder, [FreqLow FreqHigh]/(FourierSampRate/2),'bandpass');
end

%XY-Plane

%Setting the grid characteristics.
GridCharacteristicsXYList = dlmread(strcat(SimulationSpecifics2dDataDirectory,'GridCharacteristicsXYList.asc'));

DimXYx = GridCharacteristicsXYList(1,1);
DimXYy = GridCharacteristicsXYList(2,1);
CoordStepXYX = GridCharacteristicsXYList(3,1);
CoordStepXYY = GridCharacteristicsXYList(4,1);

%Filtering and collecting data.
AtmBoundListStartCoordXY = dlmread(strcat(DataStorageDirectory,'AtmBoundListStartCoordXY.asc'));
DimAtmBoundListStartCoordXY = dlmread(strcat(DataStorageDirectory,'DimAtmBoundListStartCoordXY.asc'));

iXYmin = 1;
iXYmax = DimAtmBoundListStartCoordXY;

for i=iXYmin:iXYmax
 
 x = num2str(AtmBoundListStartCoordXY(i,1));
 y = num2str(AtmBoundListStartCoordXY(i,2));
 
 %Rho
 DataListRhoSecXY = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XY-Plane/rho/DataListSecXY',x,y,'.asc'));  
 DataListRhoXY = DataListRhoSecXY(:,2);
 FilteredRhoXYData = filter(ButtB,ButtA,DataListRhoXY); 
 if FilterType == 1 
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/rhoFilt/HighFilteredDataList',x,y,'.asc'),[TimeListSec FilteredRhoXYData],'delimiter','\t');
 elseif FilterType == 2
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/rhoFilt/LowFilteredDataList',x,y,'.asc'),[TimeListSec FilteredRhoXYData],'delimiter','\t');
 elseif FilterType == 3
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/rhoFilt/BandFilteredDataList',x,y,'Peak',Peak,'.asc'),[TimeListSec FilteredRhoXYData],'delimiter','\t');
 end
 
 %Press1
 DataListPressSecXY1 = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XY-Plane/press1/DataListSecXY',x,y,'.asc'));
 DataListPressXY1 = DataListPressSecXY1(:,2);
 FilteredPressXY1Data = filter(ButtB,ButtA,DataListPressXY1); 
 if FilterType == 1 
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/press1Filt/HighFilteredDataList',x,y,'.asc'),[TimeListSec FilteredPressXY1Data],'delimiter','\t');
 elseif FilterType == 2
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/press1Filt/LowFilteredDataList',x,y,'.asc'),[TimeListSec FilteredPressXY1Data],'delimiter','\t');
 elseif FilterType == 3
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/press1Filt/BandFilteredDataList',x,y,'Peak',Peak,'.asc'),[TimeListSec FilteredPressXY1Data],'delimiter','\t');
 end
 
 %Press2
 DataListPressSecXY2 = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XY-Plane/press2/DataListSecXY',x,y,'.asc'));
 DataListPressXY2 = DataListPressSecXY2(:,2);
 FilteredPressXY2Data = filter(ButtB,ButtA,DataListPressXY2); 
 if FilterType == 1 
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/press2Filt/HighFilteredDataList',x,y,'.asc'),[TimeListSec FilteredPressXY2Data],'delimiter','\t');
 elseif FilterType == 2
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/press2Filt/LowFilteredDataList',x,y,'.asc'),[TimeListSec FilteredPressXY2Data],'delimiter','\t');
 elseif FilterType == 3
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/press2Filt/BandFilteredDataList',x,y,'Peak',Peak,'.asc'),[TimeListSec FilteredPressXY2Data],'delimiter','\t');
 end

 %Press3
 %DataListPressSecXY3 = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XY-Plane/press3/DataListSecXY',x,y,'.asc'));
 %DataListPressXY3 = DataListPressSecXY3(:,2);
 %FilteredPressXY3Data = filter(ButtB,ButtA,DataListPressXY3); 
 %if FilterType == 1 
 %dlmwrite(strcat(DataStorageDirectory,'XY-Plane/press3Filt/HighFilteredDataList',x,y,'.asc'),[TimeListSec FilteredPressXY3Data],'delimiter','\t');
 %elseif FilterType == 2
 %dlmwrite(strcat(DataStorageDirectory,'XY-Plane/press3Filt/LowFilteredDataList',x,y,'.asc'),[TimeListSec FilteredPressXY3Data],'delimiter','\t');
 %elseif FilterType == 3
 %dlmwrite(strcat(DataStorageDirectory,'XY-Plane/press3Filt/BandFilteredDataList',x,y,'Peak',Peak,'.asc'),[TimeListSec FilteredPressXY3Data],'delimiter','\t');
 %end

 %ContrVelx 
 DataListVelxSecXY = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XY-Plane/velx/DataListSecXY',x,y,'.asc'));  
 DataListVelxXY = DataListVelxSecXY(:,2);
 FilteredContrVelxXYData = filter(ButtB,ButtA,DataListVelxXY); 
 if FilterType == 1 
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/velxFilt/HighFilteredDataList',x,y,'.asc'),[TimeListSec FilteredContrVelxXYData],'delimiter','\t');
 elseif FilterType == 2
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/velxFilt/LowFilteredDataList',x,y,'.asc'),[TimeListSec FilteredContrVelxXYData],'delimiter','\t');
 elseif FilterType == 3
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/velxFilt/BandFilteredDataList',x,y,'Peak',Peak,'.asc'),[TimeListSec FilteredContrVelxXYData],'delimiter','\t');
 end
 
 %ContrVely 
 DataListVelySecXY = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XY-Plane/vely/DataListSecXY',x,y,'.asc'));  
 DataListVelyXY = DataListVelySecXY(:,2);
 FilteredContrVelyXYData = filter(ButtB,ButtA,DataListVelyXY); 
 if FilterType == 1 
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/velyFilt/HighFilteredDataList',x,y,'.asc'),[TimeListSec FilteredContrVelyXYData],'delimiter','\t');
 elseif FilterType == 2
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/velyFilt/LowFilteredDataList',x,y,'.asc'),[TimeListSec FilteredContrVelyXYData],'delimiter','\t');
 elseif FilterType == 3
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/velyFilt/BandFilteredDataList',x,y,'Peak',Peak,'.asc'),[TimeListSec FilteredContrVelyXYData],'delimiter','\t');
 end
 
 %ContrVelz 
 DataListVelzSecXY = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XY-Plane/velz/DataListSecXY',x,y,'.asc'));  
 DataListVelzXY = DataListVelzSecXY(:,2);
 FilteredContrVelzXYData = filter(ButtB,ButtA,DataListVelzXY); 
 if FilterType == 1 
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/velzFilt/HighFilteredDataList',x,y,'.asc'),[TimeListSec FilteredContrVelzXYData],'delimiter','\t');
 elseif FilterType == 2
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/velzFilt/LowFilteredDataList',x,y,'.asc'),[TimeListSec FilteredContrVelzXYData],'delimiter','\t');
 elseif FilterType == 3
 dlmwrite(strcat(DataStorageDirectory,'XY-Plane/velzFilt/BandFilteredDataList',x,y,'Peak',Peak,'.asc'),[TimeListSec FilteredContrVelzXYData],'delimiter','\t');
 end
 
end

%XZ-Plane

%Setting the grid characteristics.
GridCharacteristicsXZList = dlmread(strcat(SimulationSpecifics2dDataDirectory,'GridCharacteristicsXZList.asc'));

DimXZx = GridCharacteristicsXZList(1,1);
DimXZz = GridCharacteristicsXZList(2,1);
CoordStepXZX = GridCharacteristicsXZList(3,1);
CoordStepXZZ = GridCharacteristicsXZList(4,1);

%Filtering and collecting data.
AtmBoundListStartCoordXZ = dlmread(strcat(DataStorageDirectory,'AtmBoundListStartCoordXZ.asc'));
DimAtmBoundListStartCoordXZ = dlmread(strcat(DataStorageDirectory,'DimAtmBoundListStartCoordXZ.asc'));

iXZmin = 1;
iXZmax = DimAtmBoundListStartCoordXZ;

for i=iXZmin:iXZmax

 x = num2str(AtmBoundListStartCoordXZ(i,1));
 z = num2str(AtmBoundListStartCoordXZ(i,2));

 %Rho
 DataListRhoSecXZ = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XZ-Plane/rho/DataListSecXZ',x,z,'.asc'));  
 DataListRhoXZ = DataListRhoSecXZ(:,2);
 FilteredRhoXZData = filter(ButtB,ButtA,DataListRhoXZ); 
 if FilterType == 1 
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/rhoFilt/HighFilteredDataList',x,z,'.asc'),[TimeListSec FilteredRhoXZData],'delimiter','\t');
 elseif FilterType == 2
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/rhoFilt/LowFilteredDataList',x,z,'.asc'),[TimeListSec FilteredRhoXZData],'delimiter','\t');
 elseif FilterType == 3
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/rhoFilt/BandFilteredDataList',x,z,'Peak',Peak,'.asc'),[TimeListSec FilteredRhoXZData],'delimiter','\t');
 end
 
 %Press1
 DataListPressSecXZ1 = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XZ-Plane/press1/DataListSecXZ',x,z,'.asc'));
 DataListPressXZ1 = DataListPressSecXZ1(:,2);
 FilteredPressXZ1Data = filter(ButtB,ButtA,DataListPressXZ1); 
 if FilterType == 1 
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/press1Filt/HighFilteredDataList',x,z,'.asc'),[TimeListSec FilteredPressXZ1Data],'delimiter','\t');
 elseif FilterType == 2
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/press1Filt/LowFilteredDataList',x,z,'.asc'),[TimeListSec FilteredPressXZ1Data],'delimiter','\t');
 elseif FilterType == 3
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/press1Filt/BandFilteredDataList',x,z,'Peak',Peak,'.asc'),[TimeListSec FilteredPressXZ1Data],'delimiter','\t');
 end
 
 %Press2
 DataListPressSecXZ2 = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XZ-Plane/press2/DataListSecXZ',x,z,'.asc'));
 DataListPressXZ2 = DataListPressSecXZ2(:,2);
 FilteredPressXZ2Data = filter(ButtB,ButtA,DataListPressXZ2); 
 if FilterType == 1 
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/press2Filt/HighFilteredDataList',x,z,'.asc'),[TimeListSec FilteredPressXZ2Data],'delimiter','\t');
 elseif FilterType == 2
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/press2Filt/LowFilteredDataList',x,z,'.asc'),[TimeListSec FilteredPressXZ2Data],'delimiter','\t');
 elseif FilterType == 3
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/press2Filt/BandFilteredDataList',x,z,'Peak',Peak,'.asc'),[TimeListSec FilteredPressXZ2Data],'delimiter','\t');
 end
 
%Press3
%DataListPressSecXZ3 = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XZ-Plane/press3/DataListSecXZ',x,z,'.asc'));
%DataListPressXZ3 = DataListPressSecXZ3(:,2);
%FilteredPressXZ3Data = filter(ButtB,ButtA,DataListPressXZ3); 
%if FilterType == 1 
%dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/press3Filt/HighFilteredDataList',x,z,'.asc'),[TimeListSec FilteredPressXZ3Data],'delimiter','\t');
%elseif FilterType == 2
%dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/press3Filt/LowFilteredDataList',x,z,'.asc'),[TimeListSec FilteredPressXZ3Data],'delimiter','\t');
%elseif FilterType == 3
%dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/press3Filt/BandFilteredDataList',x,z,'Peak',Peak,'.asc'),[TimeListSec FilteredPressXZ3Data],'delimiter','\t');
%end
 
 %ContrVelx 
 DataListVelxSecXZ = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XZ-Plane/velx/DataListSecXZ',x,z,'.asc'));  
 DataListVelxXZ = DataListVelxSecXZ(:,2);
 FilteredContrVelxXZData = filter(ButtB,ButtA,DataListVelxXZ); 
 if FilterType == 1 
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/velxFilt/HighFilteredDataList',x,z,'.asc'),[TimeListSec FilteredContrVelxXZData],'delimiter','\t');
 elseif FilterType == 2
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/velxFilt/LowFilteredDataList',x,z,'.asc'),[TimeListSec FilteredContrVelxXZData],'delimiter','\t');
 elseif FilterType == 3
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/velxFilt/BandFilteredDataList',x,z,'Peak',Peak,'.asc'),[TimeListSec FilteredContrVelxXZData],'delimiter','\t');
 end 
 
 %ContrVely 
 DataListVelySecXZ = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XZ-Plane/vely/DataListSecXZ',x,z,'.asc'));  
 DataListVelyXZ = DataListVelySecXZ(:,2);
 FilteredContrVelyXZData = filter(ButtB,ButtA,DataListVelyXZ); 
 if FilterType == 1 
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/velyFilt/HighFilteredDataList',x,z,'.asc'),[TimeListSec FilteredContrVelyXZData],'delimiter','\t');
 elseif FilterType == 2
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/velyFilt/LowFilteredDataList',x,z,'.asc'),[TimeListSec FilteredContrVelyXZData],'delimiter','\t');
 elseif FilterType == 3
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/velyFilt/BandFilteredDataList',x,z,'Peak',Peak,'.asc'),[TimeListSec FilteredContrVelyXZData],'delimiter','\t');
 end 

 %ContrVelz 
 DataListVelzSecXZ = dlmread(strcat('/home/athanstav/Simulations/CygnusRuns/',SimulationName,'/',SimulationName,'Data/DataStorage/2d_DataStorage/XZ-Plane/velz/DataListSecXZ',x,z,'.asc'));  
 DataListVelzXZ = DataListVelzSecXZ(:,2);
 FilteredContrVelzXZData = filter(ButtB,ButtA,DataListVelzXZ); 
 if FilterType == 1 
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/velzFilt/HighFilteredDataList',x,z,'.asc'),[TimeListSec FilteredContrVelzXZData],'delimiter','\t');
 elseif FilterType == 2
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/velzFilt/LowFilteredDataList',x,z,'.asc'),[TimeListSec FilteredContrVelzXZData],'delimiter','\t');
 elseif FilterType == 3
 dlmwrite(strcat(DataStorageDirectory,'XZ-Plane/velzFilt/BandFilteredDataList',x,z,'Peak',Peak,'.asc'),[TimeListSec FilteredContrVelzXZData],'delimiter','\t');
 end 

end


