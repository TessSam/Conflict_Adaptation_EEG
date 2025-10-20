function place=numplace(num)
% This function is the solution for someone who has no mathmatical talent
% to come up with a one line algorithm to get the decimal place of numbers.
% When you are bored, just keep adding if statement, please.


if 0 <= num && num < 10
    place=1;%1~9
    
elseif 10 <= num && num < 100
    place =2;
    
elseif 100 <= num && num < 1000
    place =3;
    
elseif 1000 <= num && num < 10000
    place =4;
    
elseif 10000 <= num && num < 100000
    place =5;
    
else 
    place =NaN;
    
end


end