%
% Copyright 2011-2012 Ben Wojtowicz
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Affero General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Affero General Public License for more details.
%
%    You should have received a copy of the GNU Affero General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Function:    lte_calc_crc
% Description: Calculates one of the LTE CRCs.
% Inputs:      in_bits  - Input bits
%              crc_type - Which CRC to compute (8, 16,
%                         24A, or 24B)
% Outputs:     crc_bits - CRC bits
% Spec:        3GPP TS 36.212 section 5.1.1 v10.1.0
% Notes:       None
% Rev History: Ben Wojtowicz 11/18/2011 Created
%              Ben Wojtowicz 01/29/2012 Fixed license statement
%              James Peroulas 07/13/2012 Modified for use in Matlab.
function [crc_bits] = lte_calc_crc(in_bits, crc_type)
    % Check crc_type
    if(crc_type == 8)
        crc_poly = [1,1,0,0,1,1,0,1,1];
        crc_len  = 8;
    elseif(crc_type == 16)
        crc_poly = [1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1];
        crc_len  = 16;
    elseif(crc_type == '24A')
        crc_poly = [1,1,0,0,0,0,1,1,0,0,1,0,0,1,1,0,0,1,1,1,1,1,0,1,1];
        crc_len  = 24;
    elseif(crc_type == '24B')
        crc_poly = [1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1];
        crc_len  = 24;
    else
        %printf("ERROR: Invalid crc_type (%s)\n", crc_type);
        error('ERROR: Invalid crc_type');
        %crc_bits = 0;
        %return;
    end

    % Define local variables
    crc_rem   = zeros(1,crc_len+1);
    tmp_array = [in_bits, zeros(1,crc_len)];

    % Loop to calculate CRC bits
    for(n=1:length(in_bits)+crc_len)
        for(m=1:crc_len)
            crc_rem(m) = crc_rem(m+1);
        end
        crc_rem(crc_len+1) = tmp_array(n);

        if(crc_rem(1) ~= 0)
            for(m=1:crc_len+1)
                crc_rem(m) = mod(crc_rem(m)+crc_poly(m), 2);
            end
        end
    end

    for(n=1:crc_len)
        crc_bits(n) = crc_rem(n+1);
    end

