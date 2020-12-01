Init;

Factor = {AllSol AllNonSol};
Var = Data.TradeOff;
Aha = Data.AhaMoment

decalage = 1; higherBar = 1;

for f =1:length(Factor)
    
figure
VarF = Var(Factor{f});
AhaF = Aha(Factor{f});

for h = 1:2
    figure
    if h ==1
        numSub = round(length(VarF)/2);
        sub = 1:1:numSub;
    elseif h ==2
    numSub = length(VarF) - round(length(VarF)/2);
    sub = numSub:length(VarF);
    end
 cte = 0   
for s =sub(1):sub(end)
   cte = cte+1
p = numSubplots(numSub)

subplot(p(1), p(2), cte)
if f ==1
    color = 'g';
else
    color = 'r'; 
end
plot(VarF{s}, 'color', color)
 hold on
 xline(AhaF(s), 'r');
            
%     VarF_Til_Aha {sub} = VarF{sub}(1:Aha(sub))
%     VarF_afterAha{sub} = VarF{sub}(Aha(sub):length(VarF{sub}))
end

end

end

        