% Spectra Master - Version 2.0 October 2017
% Damani A. Driver, Department of Earth and Ocean Sciences,
% University of North Carolina, Wilmington:

	% The Spectra Master program uses RDSAC.m to convert .SAC time series
    % data into a MATLAB structure file 
	% storing time and signal data in 'struct'.t in 'struct'.d field 
    % vectors, respectively.
	% RDSAC.m file download and documentation are available at: 
	% https://www.mathworks.com/matlabcentral/fileexchange/46356-rdsac-and-mksac--read-and-write-sac-seismic-data-file

	% Computation of location signal average and PSD values requires that 
    % all time series data (.SAC) files be present in the same
	% working directory as RDSAC. The naming convention used to read in 
    % .SAC files is dependent on the file naming protocol used
	% in mseed2sac, which is available here:
	% https://github.com/iris-edu/mseed2sac

	% Spectral analysis is computed on a month-by-month basis
    % requiring the User to define the following variables for each 
	% month's data:
	% HHXmonth which defines the month, year, and station info
	% HHXday which defines the day and hour information
	% files which defines the specific calendar day and range 
	% where X defines the specific broadband channel being used (E N or Z)
	% these variables must be defined for each broadband channel
	
	% Further documentation can be found within program code

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Time series data are available from 10/13/2016-February 2017
% Low noise periods were estimated to fall between 3-5am

% Trial script uses October 13-31, 2016 data for demonstration purposes

%??NOTE:
%?? mseed2sac automatically appends sequential hourly series to the initial
%---  hour input
%?? eg:  mseed2sac *201610.HHE.030000.miniseed + *201610.HHE.040000.miniseed + *201610.HHE.050000.miniseed
%???? becomes *201610.HHE.030000.SAC containing 1,080,000 (360,000 x 3) data points

%----------------------------
% Read in all data by channel
%----------------------------

% For each time series, a new vector should be created containing 1,080,000 data points
% Where n = total number of time series data points
% For each month, calendar day range must be defined using variable files
% n = number of days in the range t0-tf for each broadband channel
% For this example, using October 2016 data, the range for n is  287:305 

%------------------------------
% Defining the master structure
%------------------------------

% Variable monthStruct stores matrices for each of 3 broadband channels using tempX.d structure var as input


	monthStruct=struct('HHE',zeros(1,1080000),'HHN',zeros(1,1080000),'HHZ',zeros(1,1080000));


%-HHE Input

    %?These variables must be reset for each month and channel
    % HHXmonth defines the first string for RDSAC file input
    % HHXday defines the end string for RDSAC input
    % Where X = broadband channel
    %?and cal = calendar day input string

		HHEmonth='XX.STN01..HHE.D.2016.';
		HHEday='.030000.SAC';
		files=[287:1:305];
		n=length(files);

		cal=[];
		HHEwild=[];
		HHEsac=zeros(1,1080000);	

    % Creation of RDSAC input variables
    %----------------------------------

	for i=1:n
		
		cal=num2str(files(i));

        %?Naming protocol: XX.STN01..HHE.D.2016.n.030000.SAC
        %?Iteration variable i (number) must be converted to string
        %?String value for file input must then be concatenated to contain string i in the n position
        %?This string must then be used to read in files for all vars in range n
		
        %?HHE wildcard
        %-------------
            %?mseed2sac naming protocol: 'XX.STN01..HHE.D.2016.' + 'cal' + '.030000.SAC'
            %?tempX = temporary structure storing 
		
		
		HHEwild=[HHEmonth,cal,HHEday];
		tempX=rdsac(HHEwild);
		HHEsac(i,:)=tempX.d;
		monthStruct.HHE(i,:)=HHEsac(i,:);
	end


%-HHN Input

		HHNmonth='XX.STN01..HHN.D.2016.';
		HHNday='.030000.SAC';
		cal=[];
		HHNwild=[];
		HHNsac=zeros(1,1080000);	

	for i=1:n	

		cal=num2str(files(i));
		HHNwild=[HHNmonth,cal,HHNday];
		tempX=rdsac(HHNwild);
		HHNsac(i,:)=tempX.d;
		monthStruct.HHN(i,:)=HHNsac(i,:);

	end

%-HHZ Input

		HHZmonth='XX.STN01..HHZ.D.2016.';
		HHZday='.030000.SAC';
		cal=[];
		HHZwild=[];
		HHZsac=zeros(1,1080000);	

	for i=1:n	

		cal=num2str(files(i));
		HHZwild=[HHZmonth,cal,HHZday];
		tempX=rdsac(HHZwild);
		HHZsac(i,:)=tempX.d;
		monthStruct.HHZ(i,:)=HHZsac(i,:);

    end

%------------------------
% Summation and filtering
%------------------------

% Where sumHHX = running sum vector of monthStruct.HHX matrix rows

	sumHHE=zeros(1,1080000);
	sumHHN=zeros(1,1080000);
	sumHHZ=zeros(1,1080000);
    
	HHEtemp=zeros(1,1080000);
	HHNtemp=zeros(1,1080000);
	HHZtemp=zeros(1,1080000);

%-HHE
	for i=1:n
        if i==1
            sumHHE=sumHHE+0;    
        else
            HHEtemp(i,:)=monthStruct.HHE(i,:);
            sumHHE=sumHHE(1,:)+HHEtemp(i,:);
        end
    end
    
%-HHN
	for i=1:n
        if i==1
            sumHHN=sumHHN+0;    
        else
            HHNtemp(i,:)=monthStruct.HHN(i,:);
            sumHHN=sumHHN(1,:)+HHNtemp(i,:);
        end
    end
    
%-HHZ	
    for i=1:n
        if i==1
            sumHHZ=sumHHZ+0;    
        else
            HHZtemp(i,:)=monthStruct.HHZ(i,:);
            sumHHZ=sumHHZ(1,:)+HHZtemp(i,:);
        end
    end
  

	avgHHE=sumHHE/n;
	avgHHN=sumHHN/n;
	avgHHZ=sumHHZ/n;

%----------------------------------------------
% Fast Fourier Transformation & Data Correction
%----------------------------------------------


    % Absolute Value Correction
    %--------------------------
    
    % Where N = correction factor for single-sidded fft plot

	N=ceil(length(avgHHN)/2);
	fftHHE=fft(avgHHE);
	fftHHN=fft(avgHHN);
	fftHHZ=fft(avgHHZ);
	fft_absHHE=abs(fftHHE).^2;
	fft_absHHN=abs(fftHHN).^2;
	fft_absHHZ=abs(fftHHZ).^2;

	% Frequency bin correction
    %-------------------------
    
    % Bin correction is required due to offset in fft output. MATLAB sets t0=1 while the frequency bin correction sets t0=0
	% Where 'fax' = frequency axis

	fax_bins=[0:N-1];


%-----------------------
% Power Spectral Density
%-----------------------

    % Variable Definition
    %--------------------
    
    % Where fs defines the sampling frequency in Hz
    
        fs=100; 

    % PSD Models
    %-----------
    
    % Where f = frequency bin vector generated from Welch's PSD estimate
    
        [psdHHE,psdfHHE]=pwelch(avgHHE,[],[],[],fs,'onesided');
        [psdHHN,psdfHHN]=pwelch(avgHHN,[],[],[],fs,'onesided');
        [psdHHZ,psdfHHZ]=pwelch(avgHHZ,[],[],[],fs,'onesided');

    % Min and Max Noise Models
    %-------------------------
    
        % Low-Noise Models
        %-----------------
        
        [lnmHHE,lnmfHHE]=pwelch(avgHHE,[],[],[],fs,'onesided','minhold');
        [lnmHHN,lnmfHHN]=pwelch(avgHHN,[],[],[],fs,'onesided','minhold');
        [lnmHHZ,lnmfHHZ]=pwelch(avgHHZ,[],[],[],fs,'onesided','minhold');
        
        % High-Noise Models
        %------------------
        
        [hnmHHE,hnmfHHE]=pwelch(avgHHE,[],[],[],fs,'onesided','maxhold');
        [hnmHHN,hnmfHHN]=pwelch(avgHHN,[],[],[],fs,'onesided','maxhold');
        [hnmHHZ,hnmfHHZ]=pwelch(avgHHZ,[],[],[],fs,'onesided','maxhold');
        
%-------------
% Chart Output
%-------------

    % Fourier Transformation
    %-----------------------
    
        %-HHE FFT

        subplot(3,3,1)
        plot(fax_bins(1:N),fft_absHHE(1:N))
    
        %-HHN FFT
    
        subplot(3,3,2)
        plot(fax_bins(1:N),fft_absHHN(1:N))
        
        %-HHZ FFT
    
        subplot(3,3,3)
    	plot(fax_bins(1:N),fft_absHHZ(1:N))

    % Power Spectral Density
    %-----------------------
    
        %-HHE PSD

        subplot(3,3,4)
        plot(psdfHHE,psdHHE,'k')
        hold on
        plot(lnmfHHE,lnmHHE,'m')
        
        plot(hnmfHHE,hnmHHE,'b')
        hold off
        %-HHN PSD

        subplot(3,3,5)
        plot(psdfHHN,psdHHN,'k')
        hold on
        plot(lnmfHHN,lnmHHN,'m')
        
        plot(hnmfHHN,hnmHHN,'b')
        hold off
        %-HHZ PSD

        subplot(3,3,6)
        plot(psdfHHZ,psdHHZ,'k')
        hold on
        plot(lnmfHHZ,lnmHHZ,'m')
        
        plot(hnmfHHZ,hnmHHZ,'b')
        hold off
        
    % Regression
    %----------------
    
        % FFT
        
        % PSD
        
        % HNM
        
        % LNM
        
%--------------------
% Data transformation
%--------------------

    subplot(3,3,7)
    plot(log(fax_bins(1:N)),10^100*log(fft_absHHE(1:N)))
    subplot(3,3,8)
	plot(fax_bins(1:N),fft_absHHN(1:N))
    subplot(3,3,9)
	plot(fax_bins(1:N),fft_absHHZ(1:N))

%--------------------------------------------------------------------------
% De-bugging
%--------------------------------------------------------------------------


    %-CHARTS--%
    
    % %Test HHE
    % tHHE=rdsac('XX.STN01..HHE.D.2016.287.030000.SAC');
    % dHHE=tHHE.d;
    % fHHE=abs(fft(dHHE));
    % subplot(2,3,4)
	% plot(fax_bins,fHHE)
    % %Test HHN
    % tHHN=rdsac('XX.STN01..HHN.D.2016.287.030000.SAC');
    % dHHN=tHHN.d;
    % fHHN=abs(fft(dHHN));
    % subplot(2,3,5)
	% plot(fax_bins,fHHN)
    % %Test HHZ
    % tHHZ=rdsac('XX.STN01..HHZ.D.2016.287.030000.SAC');
    % dHHZ=tHHZ.d;
    % fHHZ=abs(fft(dHHZ));
    % subplot(2,3,6)
	% plot(fax_bins,fHHZ)
