//Kenneth et Anne-Marie
//Allez dans la console scilab pour commencer
clear;

////MODELES\\\\
//1. Équation logistique
function y=F1(t,u)
    a=0.5;
    b=0.1;
    y=zeros(u);
    y=a*u*(1-b*u);
endfunction
//meme fonction mais on change les paramettres a et b
function y=F2(t,u) 
    a=7;
    b=0.14;
    y=zeros(u);
    y=a*u*(1-b*u);
endfunction
//2 Système de Lotka-Volterra
function [y]=F3(t,u)
    a=1;
    b=0.2;
    c=0.5;
    d=0.1;
    y(1)=a*u(1)-b*u(1)*u(2);
    y(2)=-c*u(2)+d*u(1)*u(2);
endfunction
//3 Modèle logistique de Verhulst
function [y]=F4(t,u)
    a=1;
    b=0.2;
    c=0.5;
    d=0.1;
    e=0.01;
    y(1)=a*u(1)-e*u(1)^2-b*u(1)*u(2);
    y(2)=-c*u(2)+d*u(1)*u(2);
endfunction
//meme fonction mais on change le paramettre e 
function [y]=F5(t,u)
    a=1;
    b=0.2;
    c=0.5;
    d=0.1;
    e=1;
    y(1)=a*u(1)-e*u(1)^2-b*u(1)*u(2);
    y(2)=-c*u(2)+d*u(1)*u(2);
endfunction
//4 Populations en compétition
function [y]=F6(t,u)
    a=1;
    th1=3/2;
    th2=1/2;
    y(1)=u(1)*(1-u(1)-th1*u(2));
    y(2)=a*u(2)*(1-u(2)-th2*u(1));
endfunction
//meme fonction mais on change les paramettres a, th1 et th2
function [y]=F7(t,u)
    a=1/2;
    th1=2;
    th2=3;
    y(1)=u(1)*(1-u(1)-th1*u(2));
    y(2)=a*u(2)*(1-u(2)-th2*u(1));
endfunction
//5 Generalisation a 3 especes 
//du modele des populations en compétition
function [y]=F8(t,u)
    a=1;
    b=2;
    //th(ij) représente l'effet de la populataion j sur la population i
    th12=3/2; 
    th21=1/2;
    th13=3/2;
    th23=1/2;
    th31=3/2;
    th32=1/2;
    y(1)=u(1)*(1-u(1)-th12*u(2)-th13*u(3));
    y(2)=a*u(2)*(1-u(2)-th21*u(1)-th23*u(3));
    y(3)=b*u(3)*(1-u(3)-th31*u(1)-th32*u(2));
endfunction

////SCHEMAS\\\\
//Schéma d’Euler explicite
function [sol]=EulerExplicite(u0,Temps,f)
    h=Temps(2)-Temps(1);
    u=u0;
    sol=u0;
    for iter=1:length(Temps)-1
        u=u+h*f(u);
        sol=[sol,u];
    end
endfunction
//Schéma d’Euler modifie
function [sol]=EulerModifie(u0,Temps,f)
    h=Temps(2)-Temps(1);
    u=u0;
    sol=u0;
    for iter=1:length(Temps)-1
        u=u+h*f(iter,u+(h/2)*f(iter,u));
        sol=[sol,u];
    end
endfunction
//RK4
function [sol]=RK4(u0,Temps,f)
    h=Temps(2)-Temps(1);
    sol=u0;
    r=u0;
    for iter=1:length(Temps)-1       
        yn1=r;
        yn2=yn1+h/2*f(t,yn1);
        yn3=yn1+h/2*f(t,yn2);
        yn4=yn1+h*f(t,yn3);
        so=yn1+h/6*(f(t,yn1)+2*f(t,yn2)+2*f(t,yn3)+f(t,yn4));
        r=so;
        sol=[sol,r]
    end
endfunction
//RK Implicite (Crank Nicholson)
function y=zer(u)
    y=u-u0-h/2*(F3(t,u0)+F3(t+h,u));
endfunction

function [sol]=CrankNicholson(u0,Temps,f)
    h=Temps(2)-Temps(1);
    sol=u0;
    r=u0;
    for iter=1:length(Temps)-1       
        yy=fsolve(r,zer);
        r=yy;
        sol=[sol,r];      
    end
endfunction

//Programme principal
//Parametres de discretisation
//Apres avoir fait les tests en temps cours,
//nous faisons passons en temps moyen (T=10)
//puis en temps long pour (T=100) pour obtenir la periodicite
T=10;
N=1000; t=linspace(0,T,N);

res=1
while(res==1)
    clc;
    clf;
    printf('=================MENU======================\n');
    printf('=1=-Équation logistique 1.(a)Euler explicite\n');
    printf('=2=-Équation logistique 1.(b)Euler modifié\n');
    printf('=3=-Système de Lotka-Volterra 2.(b)Euler explicite\n');
    printf('=4=-Système de Lotka-Volterra 2.(c)Schéma RK4\n');
    printf('=5=-Système de Lotka-Volterra 2.(d)Schéma RK implicite\n');
    printf('=6=-Modèle logistique de Verhulst 3.Euler explicite\n');
    printf('=7=-Modèle logistique de Verhulst 3.Schéma RK4\n');
    printf('=8=-Populations en compétition 4.Euler explicite\n');
    printf('=9=-Populations en compétition 4.Schéma RK4\n');
    printf('=10=-Généralisation Populations en compétition 5.Euler explicite\n');
    printf('=11=-Généralisation Populations en compétition 5.Schéma RK4\n');
    printf('-OTHER-\n');
    printf('=12=-Representation du nombre de predateur en fonction du nombre de proies(Lotka-Volterra-RK4)\n');
    printf('=13=-Representation 3D de la generalisation choisie(Competition-RK4)\n');
    choix=input('\n===========Quel est votre choix?===========\n');
    
    if(choix==1)
        //on se place avant le point d'equilibre 1/b
        //la population croit (vers le point d'equilibre)
        u0=9; //b=0.1
        subplot(2,2,1)
        solN1=EulerExplicite(u0,t,F1);
        plot(t,solN1,'r+')
        y1 = ode(u0, 0, t, F1); //solution "exacte"
        plot(t,y1)
        legend(['EulerExplicite';'solution exacte ode']);
        xtitle("Params (a=0.5 b=0.1) Val initiale (u0=9)")
        
        u0=1; //b=0.14
        subplot(2,2,2)
        solN2=EulerExplicite(u0,t,F2);
        plot(t,solN2,'r+')
        y2 = ode(u0, 0, t, F2); //solution "exacte"
        plot(t,y2)
        legend(['EulerExplicite';'solution exacte ode']);
        xtitle("Params (a=7 b=0.14) Val initiale (u0=1)")
        
        //on se place apres le point d'equilibre 1/b
        //la population decroit (vers le point d'equilibre)
        u0=11; //b=0.1
        subplot(2,2,3)
        solN1=EulerExplicite(u0,t,F1);
        plot(t,solN1,'r+')
        y1 = ode(u0, 0, t, F1); //solution "exacte"
        plot(t,y1)
        legend(['EulerExplicite';'solution exacte ode']);
        xtitle("Params (a=0.5 b=0.1) Val initiale (u0=11)")
        
        u0=50; //b=0.14
        subplot(2,2,4)
        solN2=EulerExplicite(u0,t,F2);
        plot(t,solN2,'r+')
        y2 = ode(u0, 0, t, F2); //solution "exacte"
        plot(t,y2)
        legend(['EulerExplicite';'solution exacte ode']);
        xtitle("Params (a=7 b=0.14) Val initiale (u0=50)")
        
        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end
    
    if(choix==2)
        //on se place avant le point d'equilibre 1/b
        //la population croit (vers le point d'equilibre)
        u0=1;
        subplot(2,2,1)
        solN1=EulerModifie(u0,t,F1);
        plot(t,solN1,'r+')
        y1 = ode(u0, 0, t, F1); //solution "exacte"
        plot(t,y1)
        legend(['Euler Modifie';'Solution exacte ode']);
        xtitle("Params (a=0.5 b=0.1) Val initiale (u0=1)")
        
        u0=1;
        subplot(2,2,2)
        solN2=EulerModifie(u0,t,F2);
        plot(t,solN2,'r+')
        y2 = ode(u0, 0, t, F2); //solution "exacte"
        plot(t,y2)
        legend(['Euler Modifie';'Solution exacte ode']);
        xtitle("Params (a=7 b=0.14) Val initiale (u0=1)")
        
        //on se place apres le point d'equilibre 1/b
        //la population decroit (vers le point d'equilibre)
        u0=11;
        subplot(2,2,3)
        solN1=EulerModifie(u0,t,F1);
        plot(t,solN1,'r+')
        y1 = ode(u0, 0, t, F1); //solution "exacte"
        plot(t,y1)
        legend(['Euler Modifie';'Solution exacte ode']);
        xtitle("Params (a=0.5 b=0.1) Val initiale (u0=11)")
        
        u0=50;
        subplot(2,2,4)
        solN2=EulerModifie(u0,t,F2);
        plot(t,solN2,'r+')
        y2 = ode(u0, 0, t, F2); //solution "exacte"
        plot(t,y2)
        legend(['Euler Modifie';'Solution exacte ode']);
        xtitle("Params (a=7 b=0.14) Val initiale (u0=50)")
        
        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end
    
    if(choix==3)
        //cas1
        u0=[5;2];
        subplot(2,2,1)
        solN1=EulerExplicite(u0,t,F3);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("Euler Explicite u0=[50;20]")
        subplot(2,2,2)
        y1 = ode(u0, 0, t, F3); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte u0=[50;20]")
        
        //cas2
        u0=[3;3];
        subplot(2,2,3)
        solN1=EulerExplicite(u0,t,F3);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("Euler Explicite u0=[30;30]")
        subplot(2,2,4)
        y1 = ode(u0, 0, t, F3); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte u0=[30;30]")

        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end
    
    if(choix==4)
        //cas1
        u0=[5;2];
        subplot(2,2,1)
        solN1=RK4(u0,t,F3);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("RK4 u0=[50;20]")
        subplot(2,2,2)
        y1 = ode(u0, 0, t, F3); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte u0=[50;20]")
        
        //cas2
        u0=[3;3];
        subplot(2,2,3)
        solN1=RK4(u0,t,F3);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("RK4 u0=[30;30]")
        subplot(2,2,4)
        y1 = ode(u0, 0, t, F3); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte u0=[30;30]")

        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end
    
    if(choix==5)
        //cas1
        u0=[5;2];
        subplot(2,2,1)
        solN1=CrankNicholson(u0,t,F3);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("Crank Nicholson u0=[50;20]")
        subplot(2,2,2)
        y1 = ode(u0, 0, t, F3); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte u0=[50;20]")
        
        //cas2
        u0=[3;3];
        subplot(2,2,3)
        solN1=CrankNicholson(u0,t,F3);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("Crank Nicholson u0=[30;30]")
        subplot(2,2,4)
        y1 = ode(u0, 0, t, F3); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte u0=[30;30]")

        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end
    
    if(choix==6)
        //cas1
        //e=0.01;
        //a=1;
        u0=[10;5];
        subplot(2,2,1)
        solN1=EulerExplicite(u0,t,F4);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("Euler Explicite e=0.01 u0=[100;50]")
        subplot(2,2,2)
        y1 = ode(u0, 0, t, F4); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte e=0.01 u0=[100;50]")
        
        //cas2
        //e=1;
        //a=1;
        u0=[0.1;0.5];
        subplot(2,2,3)
        solN1=EulerExplicite(u0,t,F5);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("Euler Explicite e=1 u0=[1;5]")
        subplot(2,2,4)
        y1 = ode(u0, 0, t, F5); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte e=1 u0=[1;5]")

        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end
    
    if(choix==7)
        //cas1
        //e=0.01;
        //a=1;
        u0=[10;5];
        subplot(2,2,1)
        solN1=RK4(u0,t,F4);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("RK4 e=0.01 u0=[100;50]")
        subplot(2,2,2)
        y1 = ode(u0, 0, t, F4); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte e=0.01 u0=[100;50]")
        
        //cas2
        //e=1;
        //a=1;
        u0=[0.1;0.5];
        subplot(2,2,3)
        solN1=RK4(u0,t,F5);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("RK4 e=1 u0=[1;5]")
        subplot(2,2,4)
        y1 = ode(u0, 0, t, F5); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte e=1 u0=[1;5]")
        
        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end
    
    if(choix==8)
        //cas1
        u0=[5;3];
        subplot(2,2,1)
        solN1=EulerExplicite(u0,t,F6);
        plot(t,solN1)
        legend(['Population 1';'Population 2']);
        xtitle("Euler Explicite u0=[50;30]")
        subplot(2,2,2)
        y1 = ode(u0, 0, t, F6); //solution "exacte"
        plot(t,y1)
        legend(['Population 1';'Population 2']);
        xtitle("Solution exacte u0=[50;30]")
        
        //cas2
        u0=[6;5];
        subplot(2,2,3)
        solN1=EulerExplicite(u0,t,F6);
        plot(t,solN1)
        legend(['Population 1';'Population 2']);
        xtitle("Euler Explicite u0=[60;50]")
        subplot(2,2,4)
        y1 = ode(u0, 0, t, F6); //solution "exacte"
        plot(t,y1)
        legend(['Population 1';'Population 2']);
        xtitle("Solution exacte u0=[60;50]")

        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end
    
    if(choix==9)
        //cas1
        u0=[5;3];
        subplot(2,2,1)
        solN1=RK4(u0,t,F6);
        plot(t,solN1)
        legend(['Population 1';'Population 2']);
        xtitle("RK4 u0=[50;30]")
        subplot(2,2,2)
        y1 = ode(u0, 0, t, F6); //solution "exacte"
        plot(t,y1)
        legend(['Population 1';'Population 2']);
        xtitle("Solution exacte u0=[50;30]")
        
        //cas2
        u0=[6;5];
        subplot(2,2,3)
        solN1=RK4(u0,t,F6);
        plot(t,solN1)
        legend(['Population 1';'Population 2']);
        xtitle("RK4 u0=[60;50]")
        subplot(2,2,4)
        y1 = ode(u0, 0, t, F6); //solution "exacte"
        plot(t,y1)
        legend(['Population 1';'Population 2']);
        xtitle("Solution exacte u0=[60;50]")
        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end
    
    if(choix==10)
        //cas1
        u0=[5;3;2];
        subplot(2,2,1)
        solN1=EulerExplicite(u0,t,F8);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("Euler Explicite u0=[50;30;20]")
        subplot(2,2,2)
        y1 = ode(u0, 0, t, F8); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte u0=[50;30;20]")
        
        //cas2
        u0=[5;5;5];
        subplot(2,2,3)
        solN1=EulerExplicite(u0,t,F8);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("Euler Explicite u0=[50;50;50]")
        subplot(2,2,4)
        y1 = ode(u0, 0, t, F8); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte u0=[50;50;50]")

        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end
    
    if(choix==11)
        //cas1
        u0=[5;3;2];
        subplot(2,2,1)
        solN1=RK4(u0,t,F8);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("RK4 u0=[50;30;20]")
        subplot(2,2,2)
        y1 = ode(u0, 0, t, F8); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte u0=[50;30;20]")
        
        //cas2
        u0=[5;5;5];
        subplot(2,2,3)
        solN1=RK4(u0,t,F8);
        plot(t,solN1)
        legend(['Proies';'Predateurs']);
        xtitle("RK4 u0=[50;50;50]")
        subplot(2,2,4)
        y1 = ode(u0, 0, t, F8); //solution "exacte"
        plot(t,y1)
        legend(['Proies';'Predateurs']);
        xtitle("Solution exacte u0=[50;50;50]")
        
        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end
    
    if(choix==12)
        u0=[5;2];
        solN1=RK4(u0,t,F3);
        plot(solN1(1,:),solN1(2,:))
        xtitle("Representation du nombre de predateur en fonction du nombre de proies")

        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end
    
    if(choix==13)
        u0=[5;3;2];
        solN1=RK4(u0,t,F8);        
        y1 = ode(u0, 0, t, F8); //solution "exacte"
        param3d1([solN1(1,:),y1(1,:)],[solN1(2,:),y1(2,:)],[solN1(3,:),y1(3,:)])
        xtitle("Representation 3D de la generalisation choisie(Competition-RK4) u0=[50;30;20]")

        printf('Continuer?\n');
        res=input('Entrez 1 pour continuer\nN''importe quoi d''autre pour arreter\n');
    end

end
printf('=*=*=*==*=*=*=THE=END=*=*=*==*=*=*=\n');
