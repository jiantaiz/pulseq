system = mr.opts('rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
Nx=256;
Nrep=128;

% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2,'Duration',0.1e-3, 'system', system);

% Define delays and ADC events
adc = mr.makeAdc(Nx,'Duration',3.2e-3, 'system', system,'delay',system.adcDeadTime);
delayTE=20e-3;
delayTR=1000e-3;

% Loop over repetitions and define sequence blocks
for i=1:Nrep
    seq.addBlock(rf);
    seq.addBlock(mr.makeDelay(delayTE));
    seq.addBlock(adc,mr.makeDelay(mr.calcDuration(adc)));
    seq.addBlock(mr.makeDelay(delayTR))
end

seq.write('fid.seq')       % Write to pulseq file

%% plot sequence and k-space diagrams

seq.plot('timeRange', [0 5*TR]);

% new single-function call for trajectory calculation
[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();

% plot k-spaces
time_axis=(1:(size(ktraj,2)))*sys.gradRasterTime;
figure; plot(time_axis, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points