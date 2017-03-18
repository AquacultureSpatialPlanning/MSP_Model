function out = Ei_MPAdisplace_proportional_to_fishable_patches(vector_MPA0_1fishable,Ei_original)

%this function take as inputs:
% - A vector of where you can fish (1) and not fish (0, ie MPA), and
% - A vector of the original fishing efforts. 
%...then generates output of fishing effort in fishable patches
%CHOOSE which redistribution pattern to implement
Ei_redistribute_wrt_MPAs_V=1;

%Notes: Two options below
% V1 is probably closer to optimal when there are few MPAs. It does not
% work well when high % MPAs because cramming too much effort into too few
% patches
% V2 is probably closer to optimal when there are many MPAs. It may not be
% optimal when low % MPAs because effort is not increasing in fishable
% patches near MPAs to take advantage of spillover

if Ei_redistribute_wrt_MPAs_V==1
    %V1 redistributes the displaced by MPA fishing effort to the fishable patches proportional to their effort levels.
    Esum_MPA=sum(Ei_original(vector_MPA0_1fishable==0)); %total displaced effort by MPAs
    Ei_fishable_original=Ei_original.*(vector_MPA0_1fishable==1); %generate vector that isolates just the original effort in fishable patches (set effort in MPA patches to zero) 
    Ei_fishable_original_scaled=Ei_fishable_original./sum(Ei_fishable_original); %scale Ei_fishable_original so that sums to one (use this to allocate MPA-displace effort proportionatly)
    Ei_fishable_from_MPAs=Esum_MPA.*Ei_fishable_original_scaled; %determine allocation of displaced MPA effort to fishable patches
    Ei_fishable=Ei_fishable_original+Ei_fishable_from_MPAs; %Add that allocated effort to the original effort in those fishable patches
    out=Ei_fishable; %Set the effort for the model
elseif Ei_redistribute_wrt_MPAs_V==2
    %V2 It removes the effort in the MPAs fishing effort and leave the other
    % efforts (in the non-MPAs) the same. 
    out=Ei_original.*(vector_MPA0_1fishable==1); %generate vector that isolates just the original effort in fishable patches (set effort in MPA patches to zero) 
end

