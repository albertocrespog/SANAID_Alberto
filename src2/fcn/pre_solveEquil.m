%%  -----------------------------------------------------------------------
%   Academic Stability Pro
%   Developed by Álvaro Fernández Cobo
%   On January 4th, 2016
%%   ----------------------------------------------------------------------
% This function returns the needed values of the
% incidences/incidence+position for mantenining an equilibrated flight

function [handles] = pre_solveEquil(handles)
    
    gen     = handles.model.general;    
    win     = handles.model.ala;
    hor     = handles.model.horizontal;
    can     = handles.model.canard;
    
    eq_F    = handles.auxData.eq_F;
    eq_M    = handles.auxData.eq_M;

    switch handles.model.conf

        case 'convencional_canard'
            i_val   = handles.auxData.iFix_val;
            switch handles.auxData.iFix
                case 'w'
                    syms it ic;
                    eq_Fsol =   eq_F(win.S,win.CL0,win.CLa,hor.S,hor.CL0,hor.CLa,...
                                hor.eta,0,can.S,can.CL0,can.CLa,can.eta,0,...
                                gen.Sref,gen.mtow*9.8065,gen.qinf,i_val,it,ic);
                    eq_Msol =   eq_M(win.S,win.MAC,win.Xca,win.CL0,win.CLa,win.CM0,...
                                hor.S,hor.MAC,hor.Xca,hor.CL0,hor.CLa,hor.CM0,hor.eta,...
                                0,can.S,can.MAC,can.Xca,can.CL0,can.CLa,can.CM0,...
                                can.eta,0,gen.Sref,gen.mtow*9.8065,gen.qinf,gen.Xcg,i_val,it,ic);
                    [sol1,sol2] = solve([eq_Fsol,eq_Msol],[it,ic]);
                    handles.model.horizontal.i = double(sol1);
                    handles.model.canard.i = double(sol2);
                case 'c'
                    syms it iw;
                    eq_Fsol =   eq_F(win.S,win.CL0,win.CLa,hor.S,hor.CL0,hor.CLa,...
                                hor.eta,0,can.S,can.CL0,can.CLa,can.eta,0,...
                                gen.Sref,gen.mtow*9.8065,gen.qinf,iw,it,i_val);
                    eq_Msol =   eq_M(win.S,win.MAC,win.Xca,win.CL0,win.CLa,win.CM0,...
                                hor.S,hor.MAC,hor.Xca,hor.CL0,hor.CLa,hor.CM0,hor.eta,...
                                0,can.S,can.MAC,can.Xca,can.CL0,can.CLa,can.CM0,...
                                can.eta,0,gen.Sref,gen.mtow*9.8065,gen.qinf,gen.Xcg,iw,it,i_val);
                    [sol1,sol2] = solve([eq_Fsol,eq_Msol],[it,iw]);
                    handles.model.horizontal.i = double(sol1);
                    handles.model.ala.i = double(sol2);
                case 'h'
                    syms iw ic;
                    eq_Fsol =   eq_F(win.S,win.CL0,win.CLa,hor.S,hor.CL0,hor.CLa,...
                                hor.eta,0,can.S,can.CL0,can.CLa,can.eta,0,...
                                gen.Sref,gen.mtow*9.8065,gen.qinf,iw,i_val,ic);
                    eq_Msol =   eq_M(win.S,win.MAC,win.Xca,win.CL0,win.CLa,win.CM0,...
                                hor.S,hor.MAC,hor.Xca,hor.CL0,hor.CLa,hor.CM0,hor.eta,...
                                0,can.S,can.MAC,can.Xca,can.CL0,can.CLa,can.CM0,...
                                can.eta,0,gen.Sref,gen.mtow*9.8065,gen.qinf,gen.Xcg,iw,i_val,ic);
                
                    [sol1,sol2] = solve([eq_Fsol,eq_Msol],[iw,ic]);
                    handles.model.ala.i = double(sol1);
                    handles.model.canard.i = double(sol2);
            end
   
        case 'convencional'
            % No canard
            syms iw it;
            eq_Fsol =   eq_F(win.S,win.CL0,win.CLa,hor.S,hor.CL0,hor.CLa,...
                        hor.eta,0,0,0,0,0,0,gen.Sref,gen.mtow*9.8065,...
                        gen.qinf,iw,it,0);
            eq_Msol =   eq_M(win.S,win.MAC,win.Xca,win.CL0,win.CLa,win.CM0,...
                        hor.S,hor.MAC,hor.Xca,hor.CL0,hor.CLa,hor.CM0,hor.eta,...
                        0,0,0,0,0,0,0,0,0,gen.Sref,gen.mtow*9.8065,...
                        gen.qinf,gen.Xcg,iw,it,0);
            [sol1,sol2] = solve([eq_Fsol,eq_Msol],[iw,it]);
            handles.model.ala.i = double(sol1);
            handles.model.horizontal.i = double(sol2);
        case 'canard'
            % No horizontal stabilizer
            syms iw ic;
            eq_Fsol =   eq_F(win.S,win.CL0,win.CLa,0,0,0,0,0,can.S,...
                        can.CL0,can.CLa,can.eta,0,gen.Sref,gen.mtow*9.8065,...
                        gen.qinf,iw,0,ic);
            eq_Msol =   eq_M(win.S,win.MAC,win.Xca,win.CL0,win.CLa,win.CM0,...
                        0,0,0,0,0,0,0,0,...
                        can.S,can.MAC,can.Xca,can.CL0,can.CLa,can.CM0,can.eta,0,...
                        gen.Sref,gen.mtow*9.8065,gen.qinf,gen.Xcg,iw,0,ic);
            [sol1,sol2] = solve([eq_Fsol,eq_Msol],[iw,ic]);
            handles.model.ala.i = double(sol1);
            handles.model.canard.i = double(sol2);
        case 'flyWing'
            % No canard and No horizontal stabilizer
            syms iw Xw
            eq_Fsol =   eq_F(win.S,win.CL0,win.CLa,0,0,0,0,0,0,0,0,0,0,...
                        gen.Sref,gen.mtow*9.8065,gen.qinf,iw,0,0);
            eq_Msol =   eq_M(win.S,win.MAC,Xw,win.CL0,win.CLa,win.CM0,...
                        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,gen.Sref,...
                        gen.mtow*9.8065,gen.qinf,gen.Xcg,iw,0,0);
            [sol1,sol2] = solve([eq_Fsol,eq_Msol],[iw,Xw]);
            handles.model.ala.i = double(sol1);
            handles.model.ala.Xca = double(sol2);
    end

end
