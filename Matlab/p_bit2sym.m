function psym = p_bit2sym(pbit0, pbit1, numbitpsym)
len = length(pbit0);
psym = ones(2^numbitpsym, len/numbitpsym);
switch numbitpsym
    case 1
        psym(1,:) = pbit0;
        psym(2,:) = pbit1;
    case 2
        psym(1,:) = pbit0(1:2:end).*pbit0(2:2:end);
        psym(2,:) = pbit1(1:2:end).*pbit0(2:2:end);
        psym(3,:) = pbit0(1:2:end).*pbit1(2:2:end);
        psym(4,:) = pbit1(1:2:end).*pbit1(2:2:end);
    case 4
        psym(1,:) = pbit0(1:4:end).*pbit0(2:4:end).*pbit0(3:4:end).*pbit0(4:4:end);
        psym(2,:) = pbit1(1:4:end).*pbit0(2:4:end).*pbit0(3:4:end).*pbit0(4:4:end);
        psym(3,:) = pbit0(1:4:end).*pbit1(2:4:end).*pbit0(3:4:end).*pbit0(4:4:end);
        psym(4,:) = pbit1(1:4:end).*pbit1(2:4:end).*pbit0(3:4:end).*pbit0(4:4:end);
        
        psym(5,:) = pbit0(1:4:end).*pbit0(2:4:end).*pbit1(3:4:end).*pbit0(4:4:end);
        psym(6,:) = pbit1(1:4:end).*pbit0(2:4:end).*pbit1(3:4:end).*pbit0(4:4:end);
        psym(7,:) = pbit0(1:4:end).*pbit1(2:4:end).*pbit1(3:4:end).*pbit0(4:4:end);
        psym(8,:) = pbit1(1:4:end).*pbit1(2:4:end).*pbit1(3:4:end).*pbit0(4:4:end);
        
        psym(9,:) = pbit0(1:4:end).*pbit0(2:4:end).*pbit0(3:4:end).*pbit1(4:4:end);
       psym(10,:) = pbit1(1:4:end).*pbit0(2:4:end).*pbit0(3:4:end).*pbit1(4:4:end);
       psym(11,:) = pbit0(1:4:end).*pbit1(2:4:end).*pbit0(3:4:end).*pbit1(4:4:end);
       psym(12,:) = pbit1(1:4:end).*pbit1(2:4:end).*pbit0(3:4:end).*pbit1(4:4:end);
        
       psym(13,:) = pbit0(1:4:end).*pbit0(2:4:end).*pbit1(3:4:end).*pbit1(4:4:end);
       psym(14,:) = pbit1(1:4:end).*pbit0(2:4:end).*pbit1(3:4:end).*pbit1(4:4:end);
       psym(15,:) = pbit0(1:4:end).*pbit1(2:4:end).*pbit1(3:4:end).*pbit1(4:4:end);
       psym(16,:) = pbit1(1:4:end).*pbit1(2:4:end).*pbit1(3:4:end).*pbit1(4:4:end);
   case 6
       psym(1,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
       psym(2,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
       psym(3,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
       psym(4,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
       psym(5,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
       psym(6,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
       psym(7,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
       psym(8,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
       
       psym(9,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
      psym(10,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
      psym(11,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
      psym(12,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
      psym(13,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
      psym(14,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
      psym(15,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
      psym(16,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit0(6:6:end);
      
      psym(17,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(18,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(19,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(20,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(21,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(22,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(23,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(24,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
       
      psym(25,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(26,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(27,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(28,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(29,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(30,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(31,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);
      psym(32,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit0(6:6:end);

      psym(33,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(34,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(35,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(36,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(37,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(38,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(39,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(40,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
       
      psym(41,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(42,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(43,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(44,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(45,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(46,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(47,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      psym(48,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit0(5:6:end).*pbit1(6:6:end);
      
      psym(49,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(50,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(51,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(52,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(53,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(54,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(55,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(56,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit0(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
       
      psym(57,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(58,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(59,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(60,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit0(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(61,:) = pbit0(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(62,:) = pbit1(1:6:end).*pbit0(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(63,:) = pbit0(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
      psym(64,:) = pbit1(1:6:end).*pbit1(2:6:end).*pbit1(3:6:end).*pbit1(4:6:end).*pbit1(5:6:end).*pbit1(6:6:end);
    case 8
        bitmap = de2bi(0 : ((2^numbitpsym)-1) );
        for idx = 1 : (2^numbitpsym)
            for idx1 = 1 : numbitpsym
                psym(idx,:) = psym(idx,:).*(pbit0(idx1:numbitpsym:end).*(~bitmap(idx,idx1)) + pbit1(idx1:numbitpsym:end).*bitmap(idx,idx1));
            end
        end
    otherwise
        disp('Invalid numbitpsym!!!');
end
