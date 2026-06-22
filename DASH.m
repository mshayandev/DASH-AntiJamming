%% DSH: Deterministic Synchronized Hopping for Anti-Jamming Wireless Communication
%  Monte-Carlo simulation comparing hopping strategies against three jammer models.
%
%  Strategies:  Fixed (sequential), Random, DSH (keyed deterministic PRF)
%  Jammers:     Random/Sweep, Reactive (1-slot delayed), Predictive (learning)
%  Metric:      Packet Delivery Ratio (PDR) = fraction of slots not collided.
%
%  Authors: Mohammad Shayan, Syed Ali Un Naqi Naqvi
%  NED University of Engineering and Technology, Karachi
%  Tested on MATLAB R2021b+ (and GNU Octave). No toolboxes required.
% ------------------------------------------------------------------------
clear; clc; close all; rng(42);

%% Parameters
N = 50; T = 2000; nTrials = 300; jammerBW = 1;
strategies = {'Fixed','Random','DSH'};
jammers    = {'Random','Reactive','Predictive'};
PDR = zeros(numel(strategies), numel(jammers), nTrials);

%% Monte-Carlo
for s = 1:numel(strategies)
  for j = 1:numel(jammers)
    for tr = 1:nTrials
      PDR(s,j,tr) = run_trial(strategies{s}, jammers{j}, N, T, jammerBW);
    end
  end
end
meanPDR = mean(PDR,3); stdPDR = std(PDR,0,3);

%% Console report
fprintf('\n=== Mean PDR (%%) over %d trials, N=%d, T=%d ===\n\n', nTrials, N, T);
fprintf('%-12s',''); for j=1:numel(jammers), fprintf('%-15s',jammers{j}); end; fprintf('\n');
for s=1:numel(strategies)
  fprintf('%-12s',strategies{s});
  for j=1:numel(jammers), fprintf('%6.2f +/-%4.2f  ',100*meanPDR(s,j),100*stdPDR(s,j)); end
  fprintf('\n');
end; fprintf('\n');

%% Fig 1: grouped bars
figure('Color','w','Position',[100 100 720 420]);
bar(100*meanPDR','grouped'); hold on;
ng=numel(jammers); nb=numel(strategies); gw=min(0.8,nb/(nb+1.5));
for s=1:nb
  x=(1:ng)-gw/2+(2*s-1)*gw/(2*nb);
  errorbar(x,100*meanPDR(s,:),100*stdPDR(s,:),'k','linestyle','none','HandleVisibility','off');
end
set(gca,'XTickLabel',jammers,'FontSize',11);
ylabel('Packet Delivery Ratio (%)'); xlabel('Jammer Model');
legend(strategies,'Location','southwest'); ylim([0 108]); grid on;
title('Hopping Strategy vs. Jammer Type'); saveas(gcf,'fig_pdr_bar.png');

%% Fig 2: cumulative PDR under predictive jammer (the headline result)
figure('Color','w','Position',[100 100 720 420]); hold on; cols=lines(numel(strategies));
for s=1:numel(strategies)
  c=cumulative_pdr(strategies{s},'Predictive',N,T,jammerBW);
  plot(1:T,100*c,'LineWidth',1.7,'Color',cols(s,:));
end
xlabel('Time Slot'); ylabel('Cumulative PDR (%)');
legend(strategies,'Location','east'); grid on; ylim([0 108]);
title('Cumulative PDR under Predictive (learning) Jammer'); saveas(gcf,'fig_pdr_time.png');

%% Fig 3: synchronization overhead
figure('Color','w','Position',[100 100 720 360]);
bar([0 T 0],'FaceColor',[0.85 0.33 0.31]);
set(gca,'XTickLabel',strategies,'FontSize',11);
ylabel('Control messages over T slots'); grid on;
title(sprintf('Synchronization Overhead (T=%d)',T)); saveas(gcf,'fig_sync_cost.png');

fprintf('Saved: fig_pdr_bar.png, fig_pdr_time.png, fig_sync_cost.png\n');

%% ===================== helpers =====================
function pdr = run_trial(strategy,jammer,N,T,bw)
  legit = generate_sequence(strategy,N,T);
  hits=0; obs=zeros(1,T); trans=ones(N,N);
  for t=1:T
    switch jammer
      case 'Random',    jamSet=randperm(N,bw);
      case 'Reactive',  if t==1,jamSet=randi(N);else,jamSet=obs(t-1);end
      case 'Predictive'
        if t<=2,jamSet=randi(N);else,[~,jamSet]=max(trans(obs(t-1),:));end
    end
    if any(jamSet==legit(t)),hits=hits+1;end
    obs(t)=legit(t);
    if t>1,trans(obs(t-1),legit(t))=trans(obs(t-1),legit(t))+1;end
  end
  pdr=1-hits/T;
end

function cum = cumulative_pdr(strategy,jammer,N,T,bw)
  legit=generate_sequence(strategy,N,T); obs=zeros(1,T); trans=ones(N,N); ok=0; cum=zeros(1,T);
  for t=1:T
    switch jammer
      case 'Predictive', if t<=2,jamSet=randi(N);else,[~,jamSet]=max(trans(obs(t-1),:));end
      case 'Reactive',   if t==1,jamSet=randi(N);else,jamSet=obs(t-1);end
      otherwise,         jamSet=randperm(N,bw);
    end
    if ~any(jamSet==legit(t)),ok=ok+1;end
    obs(t)=legit(t);
    if t>1,trans(obs(t-1),legit(t))=trans(obs(t-1),legit(t))+1;end
    cum(t)=ok/t;
  end
end

function seq = generate_sequence(strategy,N,T)
  seed=randi(1e9); seq=zeros(1,T);
  switch strategy
    case 'Fixed',  for t=1:T,seq(t)=mod(t-1,N)+1;end
    case 'Random', seq=randi(N,1,T);
    case 'DSH'
      K=uint32(2654435761);
      for t=1:T, seq(t)=dsh_prf(K,uint32(seed),uint32(t),N); end
  end
end

function ch = dsh_prf(K,S,t,N)
  % Overflow-safe integer PRF (xorshift-mult mixing). Deterministic given
  % (K,S,t); statistically uniform without K. Uses double() for 32-bit mults
  % to avoid uint32 saturation, then masks back to 32 bits.
  m=2^32;
  x = double(bitxor(K,S));
  x = mod(x + mod(double(t)*2246822519, m), m);
  x = mod(x*2654435761, m);
  x = bitxor(uint32(x), bitshift(uint32(x),-15));
  x = mod(double(x)*2246822519, m);
  x = bitxor(uint32(x), bitshift(uint32(x),-13));
  ch = double(mod(x, uint32(N)))+1;
end
