function MyWrite(obj,filename,Id)
%WRITE Write sequence to file.
%   WRITE(seqObj, filename) Write the sequence data to the given
%   filename using the open file format for MR sequences.
%
%   Examples:
%   Write the sequence file to the my_sequences directory
%
%       write(seqObj,'my_sequences/gre.seq')
%
% See also  read

fid=fopen(filename, 'w');
assert(fid ~= -1, 'Cannot open file: %s', filename);
fprintf(fid, '# Pulseq sequence file\n');
fprintf(fid, '# Created by MATLAB mr toolbox\n\n');

if ~isempty(obj.definitions)
    fprintf(fid, '[DEFINITIONS]\n');
    keys = obj.definitions.keys;
    values = obj.definitions.values;
    for i=1:length(keys)
        fprintf(fid, '%s ', keys{i});
        if (ischar(values{i}))
            fprintf(fid, '%s ', values{i});
        else
            fprintf(fid, '%g ', values{i});
        end
        fprintf(fid, '\n');
    end
	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, '# Format of RF events:\n');
fprintf(fid,'t Real Imag phase :\n');
%fprintf(fid, '%d %12g %d %d %g %g %g\n', [k libData1 delay libData2]);
t0=0;
for iB=1:length(obj.blockEvents)
    block=obj.getBlock(iB);
    if ~isempty(block.rf)
        rf=block.rf;
        [tc,ic]=mr.calcRfCenter(rf);
        t=rf.t+rf.delay;
        tc=tc+rf.delay;
        fprintf(fid,'%d\n',ceil(iB/4));
        amplitude=rf.signal*exp(1i*rf.phaseOffset).*exp(1i*2*pi*rf.t*rf.freqOffset);
        fprintf(fid,['%.6f,%.6f,%.6f,%.6f,%.6f\n'],transpose([t0+t*10^6,abs(amplitude),real(amplitude),imag(amplitude),angle(amplitude)]));
    end 
    if iB==Id
       amplitude=rf.signal*exp(1i*rf.phaseOffset).*exp(1i*2*pi*rf.t*rf.freqOffset);
       figure();
       subplot(411);
       plot(t0+t,real(amplitude));
       title("Real");
       subplot(412);
       plot(t0+t,imag(amplitude));
       title("Imag");
       subplot(413);
       plot(t0+t,abs(amplitude));
       title("Abs");
       subplot(414);
       plot(t0+t,angle(amplitude),...
           t0+tc,angle(rf.signal(ic)*exp(1i*rf.phaseOffset).*exp(1i*2*pi*rf.t(ic)*rf.freqOffset)),'*');
        title("Phase");
    end
    t0=t0+mr.calcDuration(block)*10^6;
    
end

% function f = plot(obj, varargin)
%             %plot Plot the sequence in a new figure.
%             %   plot(seqObj) Plot the sequence
%             %
%             %   plot(...,'Type',type) Plot the sequence with gradients
%             %   displayed according to type: 'Gradient' or 'Kspace'.
%             %
%             %   plot(...,'TimeRange',[start stop]) Plot the sequence
%             %   between the times specified by start and stop.
%             %
%             %   plot(...,'TimeDisp',unit) Display time in:
%             %   's', 'ms' or 'us'.
%             %
%             %   f=plot(...) Return the new figure handle.
%             %
%             validPlotTypes = {'Gradient','Kspace'};
%             validTimeUnits = {'s','ms','us'};
%             persistent parser
%             if isempty(parser)
%                 parser = inputParser;
%                 parser.FunctionName = 'plot';
%                 parser.addParamValue('type',validPlotTypes{1},...
%                     @(x) any(validatestring(x,validPlotTypes)));
%                 parser.addParamValue('timeRange',[0 inf],@(x)(isnumeric(x) && length(x)==2));
%                 parser.addParamValue('timeDisp',validTimeUnits{1},...
%                     @(x) any(validatestring(x,validTimeUnits)));
%             end
%             parse(parser,varargin{:});
%             opt = parser.Results;
%             
%             fig=figure;
%             if nargout>0
%                 f=fig;
%             end
%             ax=zeros(1,6);
%             for i=1:6
%                 ax(i)=subplot(3,2,i);
%             end
%             ax=ax([1 3 5 2 4 6]);   % Re-order axes
%             arrayfun(@(x)hold(x,'on'),ax);
%             arrayfun(@(x)grid(x,'on'),ax);
%             labels={'ADC','RF mag (Hz)','RF ph (rad)','Gx (kHz/m)','Gy (kHz/m)','Gz (kHz/m)'};
%             arrayfun(@(x)ylabel(ax(x),labels{x}),1:6);
%             
%             tFactorList = [1 1e3 1e6];
%             tFactor = tFactorList(strcmp(opt.timeDisp,validTimeUnits));
%             xlabel(ax(3),['t (' opt.timeDisp ')']);
%             xlabel(ax(6),['t (' opt.timeDisp ')']);
%             
%             t0=0;
%             %for iB=1:size(obj.blockEvents,1)
%             for iB=1:length(obj.blockEvents)
%                 block = obj.getBlock(iB);
%                 isValid = t0>=opt.timeRange(1) && t0<=opt.timeRange(2);
%                 
%                 if isValid
%                     if ~isempty(block.adc)
%                         adc=block.adc;
%                         t=adc.delay + (0:adc.numSamples-1)*adc.dwell;
%                         plot(tFactor*(t0+t),zeros(size(t)),'rx','Parent',ax(1));
%                     end
%                     
%                     
%                     
%                     if ~isempty(block.rf)
%                         rf=block.rf;
%                         [tc,ic]=mr.calcRfCenter(rf);
%                         t=rf.t + rf.delay;
%                         tc=tc + rf.delay;
%                         plot(tFactor*(t0+t),abs(rf.signal),'Parent',ax(2));
%                         plot(tFactor*(t0+t), angle(rf.signal    *exp(1i*rf.phaseOffset).*exp(1i*2*pi*rf.t    *rf.freqOffset)),...
%                              tFactor*(t0+tc),angle(rf.signal(ic)*exp(1i*rf.phaseOffset).*exp(1i*2*pi*rf.t(ic)*rf.freqOffset)),'xb',...
%                              'Parent',ax(3));
%                     end
%                     gradChannels={'gx','gy','gz'};
%                     for j=1:length(gradChannels)
%                         grad=block.(gradChannels{j});
%                         if ~isempty(grad)
%                             if strcmp(grad.type,'grad')
%                                 % we extend the shape by adding the first 
%                                 % and the last points in an effort of 
%                                 % making the display a bit less confusing...
%                                 t=grad.delay + [0; grad.t + (grad.t(2)-grad.t(1))/2; grad.t(end) + grad.t(2)-grad.t(1)];
%                                 waveform=1e-3* [grad.first; grad.waveform; grad.last];
%                             else
%                                 t=cumsum([0 grad.delay grad.riseTime grad.flatTime grad.fallTime]);
%                                 waveform=1e-3*grad.amplitude*[0 0 1 1 0];
%                             end
%                             plot(tFactor*(t0+t),waveform,'Parent',ax(3+j));
%                         end
%                     end                
%                 end
%                 t0=t0+mr.calcDuration(block);
%             end
%             
%             % Set axis limits and zoom properties
%             dispRange = tFactor*[opt.timeRange(1) min(opt.timeRange(2),t0)];
%             arrayfun(@(x)xlim(x,dispRange),ax);
%             linkaxes(ax(:),'x')
%             h = zoom(fig);
%             setAxesZoomMotion(h,ax(1),'horizontal');
% end


fclose(fid);
end
