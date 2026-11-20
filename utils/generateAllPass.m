% function Gap=generateAllPass(Gp)
% Gptf=tf(Gp);
% polist=pole(Gptf);
% rkpolist=find(polist>0);
% gtz=size(rkpolist);
% mz=mod(gtz(1,1),2);
% if mz==0
%     k=1;
%     %gn=1;
% else
%     k=-1;
%     %gn=zpk(1000,-1000,1);
% end
% aa=polist(rkpolist);
% Gap=zpk([aa],[-aa],k);
% %Gap=minreal(Gap*gn,1e-3);

function Gap = generateAllPass(Gp)
    Gptf = tf(Gp);
    polist = pole(Gptf);
    
    % Find unstable poles (Re > 0)
    rkpolist = find(real(polist) > 0);
    
    if isempty(rkpolist)
        % If no unstable poles, return identity system (no need for all-pass)
        Gap = tf(1);
        return;
    end
    
    % Extract unstable poles
    unstable_poles = polist(rkpolist);
    
    % Create corresponding zeros at mirrored stable locations
    stable_poles = -unstable_poles;
    
    % Define all-pass transfer function
    Gap = zpk(unstable_poles, stable_poles, 1);
end
