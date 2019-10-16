function bot = find_bottom(rcv)
% function to find the bottom depth along strike for a unicycle source obj
% Rishav Mallick, 2018, Earth Observatory of Singapore

Ustrike = unique(rcv.strike);
Udip = unique(rcv.dip);

bot = rcv.strike > 500;

if length(Ustrike) > length(Udip)
    nl = length(Ustrike);
    for i = 1:nl
        IF = rcv.strike==Ustrike(i);
        bot = bot | rcv.xc(:,3)==min(rcv.xc(IF,3));
    end        
else
    nl = length(Udip);
    for i = 1:nl
        IF = rcv.dip==Usdip(i);
        bot = bot | rcv.xc(:,3)==min(rcv.xc(IF,3));
    end 
end

end