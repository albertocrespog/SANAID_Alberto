%Este código obtiene el CYbeta de los twin vertical equivalentes para
%aplicar la formulación de Darcorporation (Help) del CY_delta_rv:

function CYbeta_twin_vertical = get_cybeta_twin(metodo,separacion)
switch separacion
    case 'no'
switch metodo
    case 'panels'
%         caso = 'T5-a0_0°-30_0m_s-Panel_TG_Twin_Vertical_equiv1.txt';
        caso = 'T5-a0_0deg-30_0m_s-Panel_TG_Twin_Vertical_equiv1.txt';
    case 'VLM2'
%         caso = 'T5-a0_0°-30_0m_s-VLM2_TG_Twin_Vertical_equiv1.txt';
        caso = 'T5-a0_0deg-30_0m_s-VLM2_TG_Twin_Vertical_equiv1.txt';
end
    case 'si'
switch metodo
    case 'panels'
%         caso = 'T5-a0_0°-30_0m_s-Panel_Twin_Vertical_equiv_sep.txt';
        caso = 'T5-a0_0deg-30_0m_s-Panel_Twin_Vertical_equiv_sep.txt';
    case 'VLM2'
%         caso = 'T5-a0_0°-30_0m_s-VLM2_Twin_Vertical_equiv_sep.txt';
        caso = 'T5-a0_0deg-30_0m_s-VLM2_Twin_Vertical_equiv_sep.txt';
end
    case 'mid'
        switch metodo
    case 'panels'
        caso = 'naara';
    case 'VLM2'
%         caso = 'T5-a0_0°-30_0m_s-VLM2_Twin_Vertical_equiv_mid.txt';
        caso = 'T5-a0_0deg-30_0m_s-VLM2_Twin_Vertical_equiv_mid.txt';
end

file_id = fopen(caso);
aux = textscan(file_id,'%n%n%n%n%n%n%n%n%n%n%n%n%n','Headerlines',8);
Data_P.beta = aux{2};
Data_P.CY    = aux{7};
fclose(file_id);
cybeta_twin = polyfit(Data_P.beta,Data_P.CY,1);
CYbeta_twin_vertical = cybeta_twin(1)*180/pi;

end