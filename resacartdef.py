# -*- coding: cp1252 -*-
from __future__ import print_function
import os
import sys
from   time  import time
import numpy as     np
import matplotlib as mpl #see: ../matplotlib/rcsetup.py
import matplotlib.pyplot as plt
from   matplotlib import cm
import matplotlib.colorbar as cb
from   matplotlib.colors import LogNorm
from   scipy.stats import norm
from   resacartparm  import *
#=====================================================================
# Traitement post param�trage
#---------------------------------------------------------------------
# D�duit du param�trage
def lidex(lin, idx) : # extrait de la list lin les �l�ments d'index idx
    lout = []
    for i in idx :
        lout.append(lin[i]);
    return lout;
if SCENARCHI>0 :
    varoutcible = lidex(varOut, targetout);
NvarIn    = len(varIn);     NvarOut= len(varOut);
IS_UVout  = False;
IS_HUVout = False;
if "SSU" in varOut and "SSV" in varOut :
    IS_UVout  = True;    
    if "SSH" in varOut :
        IS_HUVout = True;   
#
# Delta x et y des donn�es de sortie (R09 pour resac*) :
dxRout, dyRout = (ResoOut[-1]*dxR01, ResoOut[-1]*dyR01);
#
#----------------------------------------------------------------------
# Calendrier had-oc
if 1 : 
    mois  = ["oct","nov","dec","jan","feb","mar","apr","may","jun,","jul","aug","sep","oct"]; 
    nday  = [  31,   30,   31,   31,   28,   31,   30,   31,   30,    31,   31,   30,    1];
    annee = 2012;
    calendrier = [];
    for i in np.arange(len(mois)) :
        if mois[i]=="jan" :
            annee = 2013;
        for j in np.arange(nday[i]) :
            #jourj = "%d/%s/%d"%(j+1,mois[i],annee)
            jourj = "%d%s%d"%(j+1,mois[i],annee)
            calendrier.append(jourj);
    calendrier = np.array(calendrier);
#----------------------------------------------------------------------
# Calculer les Coordonn�es G�ographique en fonction des r�solutions une fois pour toute
lonfr=40.; lonto=65.; latfr=26.; latto=45.; # Coordonn�e :26�N, 45�N; 40�W, 65�W
def coordgeo(NL_,NC_,unsn) :
    lxticks = np.arange(lonto, lonfr-1, (lonfr-lonto)/(NC_-1));
    lxticks = lxticks.astype(int);      # l comme label ; arrondi au degr� le plus proche
    xxticks = np.arange(len(lxticks));  # coordonn�e en x
    #iunsn   = np.arange(0,NC_,unsn);    # prendre un point sur n - P2
    iunsn   = np.arange(0,NC_,unsn).astype(int);    # prendre un point sur n  - P3
    xxticks = xxticks[iunsn]
    lxticks = lxticks[iunsn];
    
    lyticks = np.arange(latfr, latto+1, (latto-latfr)/(NL_-1)); #NL valeurs de 25�->45� (de bas en haut)
    lyticks = lyticks.astype(int);      # l comme label ; arrondi au degr� le plus proche
    yyticks = np.arange(len(lyticks));  # coordonn�e en y
    #iunsn   = np.arange(0,NL_,unsn);    # prendre un point sur n - P2
    iunsn   = np.arange(0,NL_,unsn).astype(int);    # prendre un point sur n  - P3
    yyticks = yyticks[iunsn]
    lyticks = lyticks[iunsn];
    lxticks = -lxticks; # <- Signe moins pours les degr�s West
    return xxticks, lxticks, yyticks, lyticks;
CoordGeo = [];
ResoInterp2 = [1]; # [1] pour Interp2reso
ResoAll     = ResoIn + ResoOut + ResoInterp2;
ResoAll     = list(np.unique(ResoAll));
for i in np.arange(len(ResoAll)) :
    # D�duire NL, NC � partir de la r�solution sachant que R01 : 1296x1377
    NL_ = NLigR01 / ResoAll[i];
    NC_ = NColR01 / ResoAll[i];
    # D�terminer un nombre de points de coordonn�es
    unsn = NL_ / 8; #2      # prendre un point sur n (8 choisi en DURE)
    rc_ = coordgeo(NL_, NC_, unsn); # coordonn�es pour cette r�solution
    CoordGeo.append(rc_);
#======================================================================
# D�finitions
#----------------------------------------------------------------------
def fit01(X, gap01=0.0) :
    ''' Ram�ne les valeurs de X dans l'intervalle [0, 1] + gap01
        On retourne les valeurs (d=max(X)-min(X) et min(X/d) qui
        permettront � la fonction fit01 (ci-dessous) de
        r�aliser l'op�ration inverse.
    '''
    a = np.min(X);
    b = np.max(X)
    d = b-a; 
    Y = X / d;
    miny = np.min(Y)
    Y = Y - miny;
    Y = Y + gap01;
    coparm = ("fit01", miny, d, gap01);
    return Y, coparm
def inv_fit01(Y,coparm) :
    ''' R�alise l'op�ration inverse de la fonction fit01 ci-dessus.
    miny et d sont les param�tres qui ont �t� retourn�s par fit01.
    '''
    miny,d,gap01 = coparm[1:];
    #
    Y = Y - gap01;
    Y = Y + miny
    X = Y * d
    return X
#----------------------------------------
def refit01(X, minx0, maxx0, gap01=0.0) :
    # minx0, maxx0 : le min et le max des donn�es initiales avec
    # lesquelles le fit01 a �t� fait
    delta = maxx0 - minx0;
    Y = (X / delta) - (minx0/delta)
    Y = Y + gap01;
    return Y
def inv_refit01(Y, minx0, maxx0, gap01=0.0) :
    Y = Y - gap01;
    delta = maxx0 - minx0;
    X = (delta*Y) + minx0;
    return X
#----------------------------------------------------------------------
def codage(X,CODAGE,gap01=0.0) : 
    if CODAGE=="fit01" :
        X, coparm = fit01(X);
    else :
        raise ValueError('codage: bad code');
    return X, coparm
def decodage(Y,coparm) :
    CODAGE = coparm[0];
    if CODAGE=="fit01" :
        X = inv_fit01(Y,coparm);
    else :
        raise ValueError("decodage: code %s is unknown"%CODAGE);
    return X
def recodage(X,coparm) :
    CODAGE = coparm[0];
    if CODAGE=="fit01" :
        miny,deltax,gap01 = coparm[1:]
        Y = X / deltax;
        Y = Y - miny
        Y = Y + gap01
    else :
        raise ValueError("recodage: code %s is unknown"%CODAGE);
    return Y
#======================================================================
def showimgdata(X, Labels=None, n=1, fr=0, interp=None, cmap=CMAP_DEF, nsubl=None, 
                vmin=None, vmax=None, facecolor='w', vnorm=None, origine='lower',
                U=None, V=None, ticks=None, qscale=30) :
    if nsubl == None :
        nbsubc = np.ceil(np.sqrt(n));
        nbsubl = np.ceil(1.0*n/nbsubc);
    else :
        nbsubl = nsubl;
        nbsubc = np.ceil(1.0*n/nbsubl);
    nbsubl=int(nbsubl); #nbsubl.astype(int)
    nbsubc=int(nbsubc); #nbsubc.astype(int)
  
    if vmin is None :
        vmin = np.min(X);
    if vmax is None :
        vmax = np.max(X);
        
    N, M, P, Q = np.shape(X);      

    ISUV=False; 
    if U is not None and V is not None :
        ISUV=True;
    ISTICKS=False;
    if ticks is not None :
        ISTICKS=True;
        xxticks, lxticks, yyticks, lyticks = ticks;

    #M, P, Q = np.shape(X[0]);      
    fig, axes = plt.subplots(nrows=nbsubl, ncols=nbsubc,
                        sharex=True, sharey=True, figsize=(12,16),facecolor=facecolor)
    fig.subplots_adjust(wspace=0.1, hspace=0.3, bottom=0.0)
    ifig = 0;
    for ax in axes.flat :
        if ifig < N and ifig < n : 
            img  = X[ifig+fr]; #img = X[i];               
            if M == 1 :
                img  = img.reshape(P,Q)
            elif M != 3 :
                print("showimgdata: Invalid data dimensions image : must be : 1xPxQ or 3xPxQ")
                sys.exit(0);
            else : #=> = 3
                img  = img.transpose(1,2,0);
            if vnorm=="LogNorm" :
                ims = ax.imshow(img, norm=LogNorm(vmin=vmin, vmax=vmax),
                                interpolation=interp, cmap=cmap, origin=origine);
            else : 
                ims = ax.imshow(img, interpolation=interp, cmap=cmap, vmin=vmin,
                                vmax=vmax, origin=origine);
            if ISUV :
                u = U[ifig+fr].reshape(P, Q);
                v = V[ifig+fr].reshape(P, Q);
                ax.quiver(u, v, scale=qscale);
            if Labels is not None :
                ax.set_title(Labels[ifig+fr],fontsize=x_figtitlesize)
            if ISTICKS :
                ax.set_xticks(xxticks); ax.set_xticklabels(lxticks);
                ax.set_yticks(yyticks); ax.set_yticklabels(lyticks);
            ifig += 1;
    cbar_ax,kw = cb.make_axes([ax for ax in axes.flat],orientation="horizontal",
                             pad=0.05,aspect=30)
    fig.colorbar(ims, cax=cbar_ax, **kw);
#
#----------------------------------------------------------------------
def showdiff (Xref, Xest, relatif=False, wdif="", cmap=CMAP_DEF,
              vdmin=None, vdmax=None) :
    # les diff�rence ne sont pas (forcement) dans les m�mes echelles que les donn�es
    # Une seule image � la fois
    #print("-- in showdiff --");
    origine = ORIGINE;
    if relatif == False : 
        if wdif=="" :
            Xdif = Xref - Xest; 
            plt.imshow(Xdif, interpolation='none', vmin=vdmin, vmax=vdmax,
                       cmap=cmap, origin=origine);
        elif wdif=='log' :
            Xdif = np.log(Xref) - np.log(Xest); 
            plt.imshow(Xdif, interpolation='none', vmin=vdmin, vmax=vdmax,
                       cmap=cmap, origin=origine);
        elif wdif=='dlog' :
            Xdif = np.log(Xref - Xest); 
            plt.imshow(Xdif, interpolation='none', vmin=vdmin, vmax=vdmax, 
                       cmap=cmap, origin=origine);
        elif wdif=='dalog' :
            Xdif = np.log(np.abs(Xref - Xest)) 
            plt.imshow(Xdif, interpolation='none', vmin=vdmin, vmax=vdmax,
                       cmap=cmap, origin=origine);
    else : # Assumed relatif == True : 
        if wdif=="" :
            Xdif = (Xref - Xest) / Xref; 
            plt.imshow(Xdif, interpolation='none', vmin=vdmin, vmax=vdmax,
                       cmap=cmap, origin=origine);
        elif wdif=='log' :
            lxref = np.log(Xref);
            Xdif = (lxref - np.log(Xest)) / lxref; 
            plt.imshow(Xdif, interpolation='none', vmin=vdmin, vmax=vdmax,
                       cmap=cmap, origin=origine);
        elif wdif=='dlog' :
            Xdif = np.log((Xref - Xest)/ Xref); 
            plt.imshow(Xdif, interpolation='none', vmin=vdmin, vmax=vdmax,
                       cmap=cmap, origin=origine);
        elif wdif=='dalog' :
            Xdif = np.log(np.abs((Xref - Xest)/ Xref)); 
            plt.imshow(Xdif, interpolation='none', vmin=vdmin, vmax=vdmax,
                       cmap=cmap, origin=origine);

    # Set coord geo as ticks ...
    NL_, NC_ = np.shape(Xref);
    reso_, xxticks, lxticks, yyticks, lyticks = getrticks(NL_);
    plt.xticks(xxticks, lxticks); plt.xlabel("longitude");
    plt.yticks(yyticks, lyticks); plt.ylabel("latitude");
    titres  = "R%d(%dx%d)"%(reso_,NL_,NC_); # titre resolution
    plt.colorbar(orientation='horizontal');
    return titres;
#
#----------------------------------------------------------------------
def plotavar (Xi, titre, wnorm, wmin, wmax, cmap=CMAP_DEF, calX0i=None, Xiref=None) :
    origine = ORIGINE;
    if wnorm=="Log" :       
        wbinf=np.log(wmin); wbsup=np.log(wmax);
        plt.imshow(np.log(Xi), vmin=wbinf, vmax=wbsup, interpolation='none', cmap=cmap, origin=origine);
        titre = titre+"Log ";
    elif wnorm=="log" : # sans indiquer vmin, vmax 
        plt.imshow(np.log(Xi), interpolation='none', cmap=cmap, origin=origine);
        titre = titre+"log ";
    elif wnorm=="LogNorm" : 
        plt.imshow(Xi, norm=LogNorm(vmin=wmin, vmax=wmax), interpolation='none', cmap=cmap, origin=origine);
        titre = titre+"LogNorm ";
    else : 
        plt.imshow(Xi, interpolation='none', cmap=cmap, vmin=wmin, vmax=wmax, origin=origine);
    #
    # Set coord geo as ticks ...
    NL_, NC_ = np.shape(Xi);
    reso_, xxticks, lxticks, yyticks, lyticks = getrticks(NL_);
    plt.xticks(xxticks, lxticks); plt.xlabel("longitude");
    plt.yticks(yyticks, lyticks); plt.ylabel("latitude");
    titre = titre + "R%d(%dx%d) "%(reso_,NL_,NC_);
    #
    if calX0i is not None :
        titre = titre + calX0i;
    if Xiref is not None :
        rmsi, Nnan, inan = nanrms(Xiref, Xi);    
        titre = titre + "\nrms=%.4e "%(rmsi);    
    plt.colorbar(orientation='horizontal');
    return titre;
#--------------------------------------------------
def showsome(X, wmin=None,wmax=None,wnorm=None, cmap=CMAP_DEF, fsize=(19,8),
             calX0=None, Xref=None, wdif=None, wdifr=None, varlib="",
             suptitre="", Xtit="Est.", vdmin=None, vdmax=None) :
    #print("-- in showsome --");
    if wmin==None :
        wmin = np.min(X);
    if wmax==None :
        wmax = np.max(X);
    n = len(X);
    if SUBP_ORIENT == "horiz" :
        subl=1; subc=n;
    else :
        fsize=(8,19);
        subl=n; subc=1;
    
    plt.figure(figsize=fsize);
    for i in np.arange(n) : # VALEURS ESTIMEE
        titre = Xtit+" "+varlib+" ";
        plt.subplot(subl,subc,i+1);
        if Xref is not None :
            titre = plotavar (X[i], titre, wnorm, wmin, wmax, cmap=cmap, calX0i=calX0[i], Xiref=Xref[i]); 
        else :
            titre = plotavar (X[i], titre, wnorm, wmin, wmax, cmap=cmap, calX0i=calX0[i]); 
        plt.title(titre, fontsize=x_figtitlesize);
    plt.suptitle(suptitre,fontsize=x_figsuptitsize);

    if Xref is not None : # TRUE VALUE (valeurs de r�f�rence)
        plt.figure(figsize=fsize);
        for i in np.arange(n) :
            titre="Ref. : "+varlib+" ";
            plt.subplot(subl,subc,i+1);
            titre = plotavar (Xref[i], titre, wnorm, wmin, wmax, cmap=cmap, calX0i=calX0[i]);
            plt.title(titre, fontsize=x_figtitlesize);
        plt.suptitle(suptitre,fontsize=x_figsuptitsize);
    #
    #wdif : "": ref-est ; "log": log(ref)-log(est) ; "dlog": log(ref-est) ; "dalog" : log(abs(ref-est))
    #
    if Xref is not None : # DIFFERENCES ...
        if wdif is not None : #DIFFERENCES
            plt.figure(figsize=fsize);            
            for i in np.arange(n) :
                titre = "Dif. : "+varlib+" ";
                plt.subplot(subl,subc,i+1);                       
                titres = showdiff(Xref[i], X[i], relatif=False, wdif=wdif,
                                  cmap=cmap, vdmin=vdmin, vdmax=vdmax);
                titre  = titre + wdif + " " + titres;
                if calX0 is not None :
                    titre = titre + calX0[i];
                plt.title(titre, fontsize=x_figtitlesize);
            plt.suptitle(suptitre,fontsize=x_figsuptitsize);
        if wdifr is not None : #DIFFERENCES RELATIVES
            plt.figure(figsize=fsize);
            for i in np.arange(n) :
                titre = "Dif. Rel. : "+varlib+" ";
                plt.subplot(subl,subc,i+1);           
                titres = showdiff(Xref[i], X[i], relatif=True, wdif=wdifr,
                                  cmap=cmap, vdmin=vdmin, vdmax=vdmax);
                titre  = titre + wdifr + " " + titres;                
                if calX0 is not None :
                    titre = titre + calX0[i];
                plt.title(titre, fontsize=x_figtitlesize);
            plt.suptitle(suptitre,fontsize=x_figsuptitsize);
#----------------------------------------------------------------------
def showquivmask (Ui, Vi, qscale=None, qmask=None, qmode=None) :
    print("-- in showquivmask --");
    Nlig, Ncol = np.shape(Ui); # Une seul image a la fois, m�me shape assumed
    UU_ = np.ones((Nlig, Ncol)) * np.nan;
    VV_ = np.ones((Nlig, Ncol)) * np.nan;
    reso_    = NLigR01/Nlig; #=NColR01/Ncol
    iquv  = quvreso.index(reso_);
    if qmask is None :
        qmask  = quvmask[iquv];
    if qscale is None :
        qscale = quvscale[iquv];
    if qmode is None :
        qmode  = quvmode[iquv];
    if qmode==1 : # 'step'
        qstep = qmask;
        UU_[0:Nlig:qstep, 0:Ncol:qstep] = Ui[0:Nlig:qstep, 0:Ncol:qstep];
        VV_[0:Nlig:qstep, 0:Ncol:qstep] = Vi[0:Nlig:qstep, 0:Ncol:qstep];
    else : # moyenne' assumed
        UM_ = makemoy(Ui.reshape(1,Nlig,Ncol),qmask,qmask);
        VM_ = makemoy(Vi.reshape(1,Nlig,Ncol),qmask,qmask);
        UU_[0:Nlig:qmask , 0:Ncol:qmask] = UM_[0];
        VV_[0:Nlig:qmask , 0:Ncol:qmask] = VM_[0];
    plt.quiver(UU_,VV_, scale=qscale);    
def showquivmaskdiff (Uiref, Viref, Uiest, Viest, qscale=None,
                      qmask=None, qmode=None, relatif=False, wdif="") :
    if relatif==False :
        if wdif=="" :
            Uidif = Uiref - Uiest;
            Vidif = Viref - Viest;
        elif wdif=='log' :
            Uidif = np.log(Uiref) - np.log(Uiest);
            Vidif = np.log(Viref) - np.log(Viest);           
        elif wdif=='dlog' :
            Uidif = np.log(Uiref - Uiest);
            Vidif = np.log(Viref - Viest);
        elif wdif=='dalog' :
            Uidif = np.log(np.abs(Uiref - Uiest));
            Vidif = np.log(np.abs(Viref - Viest));
    else : # Assumed relatif == True : 
        if wdif=="" :
            Uidif = (Uiref - Uiest) / Uiref;
            Vidif = (Viref - Viest) / Uiest;
        elif wdif=='log' :
            luref = np.log(Uiref);     lvref = np.log(Viref);
            Uidif  = (luref - np.log(Uiest)) / luref;
            Vidif  = (lvref - np.log(Viest)) / lvref;           
        elif wdif=='dlog' :
            Uidif = np.log((Uiref - Uiest) / Uiref);
            Vidif = np.log((Viref - Viest) / Viref);
        elif wdif=='dalog' :
            Uidif = np.log(np.abs((Uiref - Uiest) / Uiref));
            Vidif = np.log(np.abs((Viref - Viest) / Viref));
    showquivmask(Uidif, Vidif, qscale=qscale, qmask=qmask, qmode=qmode);              
#--------------------------------------------------
def showxquiv(X, U, V, qscale=None, qmask=None, qmode=None, xvmin=None,
              xvmax=None, xvnorm=None, cmap=CMAP_DEF, fsize=(20,9), calX0=None, 
              Xref=None, Uref=None, Vref=None, wdif=None, wdifr=None, suptitre="") :
    print("-- in showxquiv --");
    # qmode : 1: 'step sinon 'moyenne'
    # qmask : Si qmode!=1, en attendant de revoir la fonction makemoy,
    #         ce doit etre un diviseur de Nlig et Ncol - not checked
    Nimg, Nlig, Ncol = np.shape(X); # same for U et V not checked
    reso_, xxticks, lxticks, yyticks, lyticks = getrticks(Nlig); # reso et coord g�o
    origine = ORIGINE;
    #
    if SUBP_ORIENT == "horiz" :
        subl=1; subc=Nimg;
    else :
        fsize=(9, 18);
        subl=Nimg; subc=1;
    #    
    # VALEURS ESTIMEES (h+quiv(uv))
    plt.figure(figsize=fsize);
    for i in np.arange(Nimg) :
        titre="Est.. ";
        plt.subplot(subl,subc,i+1);
        plt.imshow(X[i], interpolation='none', cmap=cmap, vmin=xvmin, vmax=xvmax, origin=origine);
        plt.colorbar(orientation='horizontal');
        showquivmask(U[i], V[i], qscale=qscale, qmask=qmask, qmode=qmode);        
        if 1 : # Avoir d�finit les coordonn�es g�ographique par r�solution une fois pour toute
            plt.xticks(xxticks, lxticks); plt.xlabel("longitude");
            plt.yticks(yyticks, lyticks); plt.ylabel("latitude");
            titre = titre + "R%d(%dx%d)\n"%(reso_,Nlig,Ncol);
        if calX0 is not None :
            titre = titre + calX0[i];
        if Xref is not None : # On suppose que c'est idem pour Uref et Vref
            RMSX, Nnan, inan = nanrms(Xref[i], X[i]);
            RMSU, Nnan, inan = nanrms(Uref[i], U[i]);
            RMSV, Nnan, inan = nanrms(Vref[i], V[i]);
            titre = titre + " rmsH=%.4e, \nrmsU=%.4e, rmsV=%.4e"%(RMSX,RMSU,RMSV);
        plt.title(titre, fontsize=x_figtitlesize);
    plt.suptitle(suptitre,fontsize=x_figsuptitsize);
    if Xref is not None : # On suppose que c'est idem pour Uref et Vref
        # VALEUR de REFERENCE (True) (h+quiv(uv))
        plt.figure(figsize=fsize);
        for i in np.arange(Nimg) :
            titre="Ref.. ";
            plt.subplot(subl,subc,i+1);
            plt.imshow(Xref[i], interpolation='none', cmap=cmap, vmin=xvmin, vmax=xvmax, origin=origine);
            plt.colorbar(orientation='horizontal');
            showquivmask(Uref[i], Vref[i], qscale=qscale, qmask=qmask, qmode=qmode);

            # Set coord geo as ticks ...
            plt.xticks(xxticks, lxticks); plt.xlabel("longitude");
            plt.yticks(yyticks, lyticks); plt.ylabel("latitude");
            titre = titre + "R%d(%dx%d)"%(reso_,Nlig,Ncol);
                
            if calX0 is not None :
                titre = titre + calX0[i];
            plt.title(titre, fontsize=x_figtitlesize);
        plt.suptitle(suptitre,fontsize=x_figsuptitsize);
        #
        if wdif is not None :   # DIFFERENCES (h+quiv(uv))
            ivar_ = tvwmm.index("SSH");
            vdmin = wdmin[ivar_];
            vdmax = wdmax[ivar_];
            plt.figure(figsize=fsize);
            for i in np.arange(Nimg) :
                titre="Dif.. ";
                plt.subplot(subl,subc,i+1);
                titres = showdiff(Xref[i], X[i], cmap=cmap, relatif=False,
                                  wdif=wdif, vdmin=vdmin, vdmax=vdmax);
                showquivmaskdiff (Uref[i], Vref[i], U[i], V[i],
                          qscale=qscale, qmask=qmask, qmode=qmode, relatif=False, wdif=wdif);
                titre = titre + " " + titres;
                if calX0 is not None :
                    titre = titre + calX0[i];
                plt.title(titre, fontsize=x_figtitlesize);
            plt.suptitle(suptitre,fontsize=x_figsuptitsize);

        if wdifr is not None :  # DIFFERENCES RELATIVES (h+quiv(uv))
            plt.figure(figsize=fsize);
            for i in np.arange(Nimg) :
                titre="Dif.. Rel. SSH +(U,V) ";
                plt.subplot(subl,subc,i+1);
                titres = showdiff(Xref[i], X[i], cmap=cmap, relatif=True,
                                  wdif=wdifr, vdmin=vdmin, vdmax=vdmax);
                showquivmaskdiff (Uref[i], Vref[i], U[i], V[i], qscale=qscale,
                                  qmask=qmask, qmode=qmode, relatif=True, wdif=wdifr);                                        
                titre = titre + " " + titres;
                if calX0 is not None :
                    titre = titre + calX0[i];
                plt.title(titre, fontsize=x_figtitlesize);
            plt.suptitle(suptitre,fontsize=x_figsuptitsize);
#--------------------------------------------------
def showhquiv(Voutb, im2show, qscale=None, qmask=None, qmode=None,
              fsize=(20,9), calX0=None, Yref=None, wdif=None, wdifr=None, suptitre="") : 
    #print("-- in showhquiv --");
    if calX0 is not None :
        calX0_ = calX0[im2show];
    iu_ = varOut.index("SSU");
    iv_ = varOut.index("SSV");
    ih_ = varOut.index("SSH",np.min([iu_,iv_])-1); # le dernier h avant u et v        
    H_  = Voutb[ih_][im2show,0,:,:];
    U_  = Voutb[iu_][im2show,0,:,:];
    V_  = Voutb[iv_][im2show,0,:,:];
    wk_ = tvwmm.index(varOut[ih_]);
    Href_=None; Uref_=None; Vref_=None;
    if Yref is not None :
        Href_ = Yref[ih_][im2show,0,:,:];
        Uref_ = Yref[iu_][im2show,0,:,:];
        Vref_ = Yref[iv_][im2show,0,:,:];        
    showxquiv(H_, U_, V_, qscale=qscale, qmask=qmask, qmode=qmode, 
              xvmin=wbmin[wk_], xvmax=wbmax[wk_], xvnorm=wnorm[wk_], cmap=wcmap[wk_],
              fsize=fsize, calX0=calX0_, Xref=Href_, Uref=Uref_, Vref=Vref_,
              wdif=wdif, wdifr=wdifr, suptitre=suptitre);
    return ih_;
#----------------------------------------------------------------------
def rms(Xref, Xest) :
    '''Calcule et retourne la RMS entre Xref et Xest
    '''
    return np.sqrt(np.sum((Xref-Xest)**2) / np.prod(np.shape(Xref)));
def rmsrel(Xref, Xest) :
    '''Calcule et retourne la RMS relative entre Xref et Xest '''
    Nall   = np.prod(np.shape(Xref));
    Errrel = (Xref - Xest) / Xref;   
    Errrel = np.sqrt(np.sum(Errrel**2)/Nall); 
    return Errrel;
def nanrms(X, Y) :
    Nitem_ = np.prod(np.shape(X))
    X_     = np.reshape(X, Nitem_);
    Y_     = np.reshape(Y, Nitem_);
    ixnan  = np.where(np.isnan(X_))[0];
    iynan  = np.where(np.isnan(Y_))[0];
    inan   = np.union1d(ixnan, iynan);
    Nnan   = len(inan)
    Nitem_ = Nitem_- Nnan;       
    RMS    = np.sqrt(np.nansum((X_-Y_)**2) / Nitem_ );
    return RMS, Nnan, inan
#----------------------------------------------------------------------
def makemoy(XB, ml=3, mc=3) :
    N,nl,nc = np.shape(XB);
    xm      = np.zeros((N,nl//ml,nc//mc));
    ii = 0;
    for i in np.arange(0,nl,ml) :
        jj = 0;
        for j in np.arange(0,nc,mc) :
            if 0 :
                moy = np.mean(XB[:,i:i+ml,j:j+mc]);
                xm[0,ii,jj] = moy; 
            else :              
                moy = np.mean(XB[:,i:i+ml,j:j+mc], axis=(1,2));
                xm[:,ii,jj] = moy;
            jj=jj+1;
        ii=ii+1;
    return xm
#-------------------------------------------------------------
def isetalea (Nimg, pcentSet) :
    pcentA, pcentV, pcentT = pcentSet;
    Ialea = np.arange(Nimg);
    np.random.seed(0); # Pour reproduire la m�me chose � chaque fois
    np.random.shuffle(Ialea);
    indA = Ialea[0:int(Nimg*pcentA)];
    indV = Ialea[int(Nimg*pcentA): int(Nimg*(pcentA+pcentV))];
    indT = Ialea[int(Nimg*(pcentA+pcentV)): int(Nimg*(pcentA+pcentV+pcentT))];
    return indA, indV, indT;    
#----------------------------------------
def splitset (Vin_brute, Vout_brute, pcentSet) :
    Nimg_ = len(Vin_brute[0]); # == len(Vou_brute[0]); not checked
    indA, indV, indT = isetalea(Nimg_, pcentSet);
    #
    VAout_brute = []; VVout_brute = []; VTout_brute = [];
    for i in np.arange(len(Vout_brute)) : # Pour chaque variable (i.e. liste)
        VAout_brute.append(Vout_brute[i][tuple([indA])]);               
        VVout_brute.append(Vout_brute[i][tuple([indV])]);               
        VTout_brute.append(Vout_brute[i][tuple([indT])]);               
    VAin_brute = []; VVin_brute = []; VTin_brute = [];
    for i in np.arange(len(Vin_brute)) : # Pour chaque variable (i.e. liste)
        VAin_brute.append(Vin_brute[i][tuple([indA])]);               
        VVin_brute.append(Vin_brute[i][tuple([indV])]);               
        VTin_brute.append(Vin_brute[i][tuple([indT])]);
    return VAin_brute, VAout_brute, VVin_brute, VVout_brute, VTin_brute, VTout_brute;
#-------------------------------------------------------------
def setresolution(VA_brute,VV_brute,VT_brute,varlue,ResoIn,ResoOut) :           
    # Make resolution for IN and OUT
    print("... making V*out_Brute");
    VAout_brute = []; VVout_brute = []; VTout_brute = [];
    for i in np.arange(NvarOut) : #varOut ['SSH', 'SSH', 'SSU', 'SSV']
        idvar = varlue.index(varOut[i])
        dvar  = VA_brute[idvar];
        if ResoOut[i] > 1 :
            dvar  = makemoy(dvar, ResoOut[i], ResoOut[i]);
        VAout_brute.append(dvar);            
        dvar  = VV_brute[idvar];
        if ResoOut[i] > 1 :
            dvar  = makemoy(dvar, ResoOut[i], ResoOut[i]);
        VVout_brute.append(dvar);            
        dvar  = VT_brute[idvar];
        if ResoOut[i] > 1 :
            dvar  = makemoy(dvar, ResoOut[i], ResoOut[i]);
        VTout_brute.append(dvar);            
    print("... making V*in_Brute");
    VAin_brute = []; VVin_brute = []; VTin_brute = [];
    for i in np.arange(NvarIn) : #varIn['SSH', 'SST', 'SST']
        idvar = varlue.index(varIn[i])
        dvar  = VA_brute[idvar];
        if ResoIn[i] > 1 :
            dvar  = makemoy(dvar, ResoIn[i], ResoIn[i]);
        VAin_brute.append(dvar);            
        dvar  = VV_brute[idvar];
        if ResoIn[i] > 1 :
            dvar  = makemoy(dvar, ResoIn[i], ResoIn[i]);
        VVin_brute.append(dvar);            
        dvar  = VT_brute[idvar];
        if ResoIn[i] > 1 :
            dvar  = makemoy(dvar, ResoIn[i], ResoIn[i]);
        VTin_brute.append(dvar);
    return VAout_brute, VVout_brute, VTout_brute, VAin_brute, VVin_brute, VTin_brute;
#-------------------------------------------------------------
def statibase(Xi, fmt=None) : # stat de base
    if fmt is None : 
        print("min=%.4f, max=%.4f, moy=%.4f, std=%.4f"
                  %(np.min(Xi), np.max(Xi), np.mean(Xi), np.std(Xi)));            
    else : 
        print("min=%.4e, max=%.4e, moy=%.4e, std=%.4e"
                  %(np.min(Xi), np.max(Xi), np.mean(Xi), np.std(Xi)));            
def stat2base (X, nmvar) : 
    for i in np.arange(len(nmvar)) :           
        print("%s : "%nmvar[i], end=''); statibase(X[i])
#-------------------------------------------------------------
def linr2 (x, y) :
    ''' b0,b1,s,R2,sigb0,sigb1 = .linreg(x,y)
    | Calcule la r�gression lin�aire de x par rapport a y.
    | En sortie :
    | b0, b1 : Coefficients de la droite de regression lineaire : y=b0+b1*x
    | R2     : Coefficient de d�termination
    '''
    N       = x.size;
    xmean   = np.mean(x);
    xmean2  = xmean*xmean;
    xcentre = x-xmean;
    ymean   = np.mean(y);
    ycentre = y-ymean;
    b1 = np.sum(ycentre*xcentre) / (np.sum(x*x) - N*xmean2);
    b0 = ymean - b1*xmean;  
    yc = b0 + b1*x;
    R2 = np.sum(pow((yc-ymean),2))/np.sum(pow(ycentre,2));  
    return b0,b1,R2
#--------------------------------------------------
def histodbl(Xref, Xest, legende, nbins=NBINS) :
    # Histogrammes double (supperpos�s) pour comparaison 
    mmin = np.min([Xref, Xest]); 
    mmax = np.max([Xref, Xest]); 
    epsi  = abs(mmax)/100000;
    mmaxX = mmax + epsi; 
    bins  = np.arange(mmin, mmaxX, (mmax-mmin)/nbins);
    bins[-1] = bins[-1] + epsi;
    plt.figure();
    plt.hist(Xref, bins=bins);
    plt.hist(Xest, bins=bins, color=[0.60,1.00,0.60,0.7]);        
    plt.legend(legende, framealpha=0.5);
    return
#--------------------------------------------------
def scatplot (Xref, Xest, ident=True, regr2=True, mksz=3.0, linewidth=2.0) :
    # Scatter plot
    mpl.rcParams['agg.path.chunksize'] = 1000000;
    plt.figure();
    plt.plot(Xref, Xest, '.g', markersize=mksz); 
    minim=np.nanmin([Xref, Xest]); maxim=np.nanmax([Xref, Xest]);
    plt.xlabel("NATL60"); plt.ylabel("RESAC");
    if ident : # ligne de l'identit�
        plt.plot([minim, maxim], [minim, maxim], '-m',linewidth=linewidth);
    if regr2 : # r�gression + R2
        b0, b1, R2 = linr2(Xref, Xest);
        xmin=min(Xref); xmax=max(Xref);
        xr = np.arange(xmin, xmax,(xmax-xmin)/200); #[xmin:(xmax-xmin)/200:xmax]';  
        yr = b1*xr + b0;
        plt.plot(xr, yr, '-b',linewidth=2.0);
        plt.legend(["", "id", "rl"], numpoints=1, loc='upper left',framealpha=0.5);
        return b0, b1, R2;
    return
#--------------------------------------------------
def setbins(mmax, mmin, nbins) :
    epsi  = abs(mmax)/100000;
    mmaxX = mmax + epsi; 
    bins  = np.arange(mmin, mmaxX, (mmax-mmin)/nbins);
    bins[-1] = bins[-1] + epsi;
    return bins
def getrticks (Nlig) :
    # get r�solution et ticks
    reso    = NLigR01/Nlig; #=NColR01/NC_
    icoord  = ResoAll.index(reso);
    xxticks, lxticks, yyticks, lyticks = CoordGeo[icoord];
    return reso, xxticks, lxticks, yyticks, lyticks    
def nl2limlat (Nlig, limlat) :
    # Nombre de ligne jusqu'� une latitude limite
    #lonfr=40.; lonto=65.;latfr=26.; latto=45.; # Coordonn�e :26�N, 45�N; 40�W, 65�W
    Nlat   = latto-latfr+1;              # =45-26+1    = 20.0
    Nllat  = Nlig / Nlat;                # =144 / Nlat = 7.20 : nombre de ligne par lattitude
                                         # v�rif : (45 - 26 +1) * 7.20 = 144
    #limlat = 35;                        # La limite de latitude choisie EN DURE
    NL_lim = int((limlat-latfr+1)*Nllat);# =(35-26+1)*7.20 = 72 : nombre de lignes jusqu'� la lattitude limite
    return NL_lim
def splitns (X, NL_lim) :
    # Split Nord-Sud
    if len(np.shape(X))==3 : # Cas Enstrophie
        pipo, Nlig, pipo = np.shape(X);
        XS = X[:,0:NL_lim,:];
        XN = X[:,NL_lim:Nlig,:];
    else : # assume==4 (cas Energie)
        pipo, pipo, Nlig, pipo = np.shape(X);
        XS = X[:,:,0:NL_lim,:];
        XN = X[:,:,NL_lim:Nlig,:];
    return XN, XS;
#--------------------------------------------------
def result(lib, strset, varname, Yref, Yest, flag_rms=0, flag_scat=0,
           flag_histo=0, nbins=NBINS, calX=None) :
    print("-- in result --"); #print(np.shape(Yref)); #(55L, 1L, 144L, 153L)
    # Au moins toujours la RMS qui est affich�e par ailleurs
    N_ = len(Yref); Allrmsi = [];
    for i in np.arange(N_) :
        rmsi, Nnan, inan = nanrms(Yref[i], Yest[i]);
        Allrmsi.append(rmsi);
    moyAllrmsi = np.mean(Allrmsi);
    titres_ = "%s %s %s, RMS by image, Mean : %.4f\n min=%.4f, max=%.4f, std=%.4f" \
               %(strset, varname, lib, moyAllrmsi, np.min(Allrmsi), np.max(Allrmsi), np.std(Allrmsi));
    nom_save="%s_%s_%s_RMS_BY_IMAGE"%(strset, varname, lib)
    print(titres_);
    #
    if flag_rms and FIGBYIMGS : # figure RMS par image
        if calX==None : # plot non tri� sur la date
            plt.figure(); plt.plot(Allrmsi, '-*');
            plt.xticks(np.arange(N_), calX[0], rotation=35, ha='right');
        else : # plot tri� sur la date
            ids = np.argsort(calX[1]); # indice de tri des dates              
            plt.figure();
            plt.plot(np.array(Allrmsi)[ids], '-*');
            plt.xticks(np.arange(N_),calX[0][ids], rotation=35, ha='right');
            #horizontalalignment='right', verticalalignment='baseline')          
        plt.plot([0,N_-1],[moyAllrmsi, moyAllrmsi]); # Trait de la moyenne
        plt.title(titres_, fontsize=x_figtitlesize);
        plt.savefig("Images/"+nom_save+".png")       
    #
    Nall_ = np.prod(np.shape(Yest));
    Yref_ = Yref.ravel();
    Yest_ = Yest.ravel();
    #                  
    if flag_scat :     # Scatter plot
        pipo, pipo, R2 = scatplot (Yref_, Yest_);
        plt.title("%s %s %s scatter plot (%d pixels)\n mean(rms)=%.6f, R2=%f"
                  %(strset, varname, lib, Nall_, moyAllrmsi, R2), fontsize=x_figtitlesize);
        nom_save2="%s_%s_%sscatter_plot_(%d_pixels)"%(strset, varname, lib, Nall_)
        plt.savefig("Images/"+nom_save2+"scatt.png")

    if flag_histo : # Histogramme des diff�rences
        Ydif_ = Yref_-Yest_;
        wk_   = tvwmm.index(varname);
        # Positionnement des min et max d'erreur (Forcage des bins sur EXP1)
        mmin = whdmin[wk_];  mmax = whdmax[wk_]; 
        bins = setbins(mmax, mmin, nbins);
        
        if 0 : # Histo non Normalis�
            plt.figure();
            Nbmax_= whnmax[wk_];
            BI = plt.hist(Ydif_, bins=bins); #print("max(BI[0])=%d"%(np.max(BI[0])))
            plt.axis([mmin, mmax, 0, Nbmax_]);
            plt.title("%s %s %s Histo des differences (Ref.-Est.)\n(%dbins)(%d pixels(%dnans) %dloss), mean(rms)=%.6f"
                      %(strset, varname, lib, nbins, Nall_, Nnan, Nall_-sum(BI[0]), moyAllrmsi), fontsize=x_figtitlesize);
            nom_save3="%s_%s_%s_Histo_des_differences"%(strset, varname, lib)
            plt.savefig("Images/"+nom_save3+"histo_non_norm.png")

        # Histogramme Normalis�
        plt.figure();
        #strset may be : "APP", "TEST", "APP-Nord" ou "APP-Sud" ou "TEST-Nord" ou "TEST-Sud"
        if len(strset)>4 :          # Pour savoir si on est dans un cas Nord-Sud ou pas, on
            NBmax_ = whSmax[wk_];   # utilise un test sur la longueur de strset ; ok c'est pas
        else :                      # top ; p'tet qu'un jour on c�era un param�tre d�di� ...
            NBmax_ = whNmax[wk_];  
        BI = plt.hist(Ydif_, bins=bins,density=True)#, normed=True);
        dbin   = bins[1] - bins[0]; 
        dhbin  = dbin / 2; 
        Centre = bins[0:nbins]+dhbin;
        ymoy   = np.mean(Ydif_);
        ystd   = np.std(Ydif_);
        y      = norm.pdf(Centre, ymoy, ystd);
        plt.plot(Centre, y, '-*r', linewidth=3);
        plt.axis([mmin, mmax, 0, NBmax_]);
        plt.title("%s %s %s Histo Normalise des Diff.(Ref.-Est.)\n mean=%f, std=%f, mean(rms)=%.6f"
                  %(strset, varname, lib, np.mean(Ydif_), np.std(Ydif_), moyAllrmsi));
        nom_save4="%s_%s_%s_Histo_Normalise_des_Diff"%(strset, varname, lib)
        plt.savefig("Images/"+nom_save4+"histo_norm.png")                    
#
#----------------------------------------------------------------------      
def distcosine2D (Wref, West) : # Matrix way
    WW    = np.sum(Wref*West, axis=1);   # les produits scalaires <wref, west>
    Nwref = np.sqrt(np.sum(Wref**2, 1)); # Normes des wref
    Nwest = np.sqrt(np.sum(West**2, 1)); # Normes des west
    costheta  = (WW/(Nwref*Nwest));
    moycosine = np.mean(costheta);
    return 1 - moycosine;
def dcosine2D(Uref,Vref,Uest,Vest) :
    Nall = np.prod(np.shape(Uref)); #idem for others not cheked
    Wref = np.array([Uref.reshape(Nall), Vref.reshape(Nall)]);
    West = np.array([Uest.reshape(Nall), Vest.reshape(Nall)]);
    return distcosine2D(Wref, West);
#----------------------------------------------------------------------
def enstrophie2d (U, V, dx, dy) :
    # U, V de la forme (Nimg, Nlig, Ncol) (16L, 144L, 153L)
    dedx = dx*2; dedy = dy*2;
    lenshape = len(np.shape(U));
    if lenshape == 3 :
        Nimg, Nlig, Ncol = np.shape(U); # same for V; not checked
        U_= U; V_= V;
    elif lenshape == 4 :
        Nimg, Ncan, Nlig, Ncol = np.shape(U); # same for V; not checked
        U_ = U.reshape(Nimg, Nlig, Ncol);
        V_ = V.reshape(Nimg, Nlig, Ncol);
        
    axX = np.arange(Ncol-2)+1;
    axY = np.arange(Nlig-2)+1;
    Enst= np.zeros((Nimg, Nlig, Ncol));
    for l in axY : # lig
        for c in axX : # col          
             dvdxdudy    = ((V_[:,l,c+1] - V_[:,l,c-1]) / dedx) \
                         - ((U_[:,l-1,c] - U_[:,l+1,c]) / dedy);
             Enst[:,l,c] = (dvdxdudy**2) / 2;             
    Enst = Enst[:, 1:Nlig-1, 1:Ncol-1]; 
    return Enst;
def divergence2d (U, V, dx, dy) :
    # U, V de la forme (Nimg, Nlig, Ncol) (16L, 144L, 153L)
    dedx = dx*2; dedy = dy*2;
    lenshape = len(np.shape(U));
    if lenshape == 3 :
        Nimg, Nlig, Ncol = np.shape(U); # same for V; not checked
        U_= U; V_= V;
    elif lenshape == 4 :
        Nimg, Ncan, Nlig, Ncol = np.shape(U); # same for V; not checked
        U_ = U.reshape(Nimg, Nlig, Ncol);
        V_ = V.reshape(Nimg, Nlig, Ncol);
        
    axX = np.arange(Ncol-2)+1;
    axY = np.arange(Nlig-2)+1;
    Div = np.zeros((Nimg, Nlig, Ncol));
    for l in axY : # lig
        for c in axX : # col            
             dudxdvdy    = ((U_[:,l,c+1] - U_[:,l,c-1]) / dedx) \
                         + ((V_[:,l-1,c] - V_[:,l+1,c]) / dedy);            
             Div[:,l,c] = (dudxdvdy**2) / 2;             
    Div = Div[:, 1:Nlig-1, 1:Ncol-1];
    return Div;
#----------------------------------------------------------------------
def phistuff (phy_ref, phy_est, varname, varlib, strset, calX, wdif,
              flag_histo, flag_scat, nbins=NBINS) :
    print("-- in phistuff --");       
    if len(np.shape(phy_ref))==4 : # cas energie
         Nimg_, Ncan_, Nlig_, Ncol_ = np.shape(phy_ref);
         phi_ref = phy_ref.reshape(Nimg_, Nlig_, Ncol_); 
         phi_est = phy_est.reshape(Nimg_, Nlig_, Ncol_);
    else : # assumed==3 : cas enstrophie
         Nimg_, Nlig_, Ncol_ = np.shape(phy_ref);
         phi_ref = phy_ref 
         phi_est = phy_est    
    Nall_ = Nimg_*Nlig_*Ncol_;   
    reso_ = NLigR01/Nlig_;
    wk_   = tvwmm.index(varname);
    #
    if strset[0]=='A' : #"APP" :
        im2show = im2showA;
    elif strset[0]=='T'  : #"TEST" :
        im2show = im2showT;
    else :
        raise ValueError("salperlipopette");
    #
    if not FIGBYPASS : 
        suptitre = "Some %s %s %s"%(strset,varlib,im2show);
        ivar_ = tvwmm.index(varname);
        vdmin = wdmin[ivar_];
        vdmax = wdmax[ivar_];
        if len(im2show) and reso_ in ResoAll :
            # ps : dans le cas Nord-Sud reso_ ne sera pas dans ResoAll
            #      et on ne devrait pas avoir besoin des images dans ce cas
            showsome(phi_est[im2show], wmin=wbmin[wk_], wmax=wbmax[wk_],
                     wnorm=wnorm[wk_], cmap=wcmap[wk_], calX0=calX[0][im2show], Xref=phi_ref[im2show],
                     wdif=wdif, varlib=varlib, suptitre=suptitre, vdmin=vdmin, vdmax=vdmax);               

    # Les Erreurs Relatives : toujours car la moyenne est utilis�e apr�s 
    RelErr_ = []; RMRelErr_ = [];                    
    for j in np.arange(Nimg_) :
        Ej_ref   = np.sum(phi_ref[j]);  
        Ej_est   = np.sum(phi_est[j])
        RelErrj_ = np.abs(Ej_ref - Ej_est) / Ej_ref;
        RelErr_.append(RelErrj_);
    MoyRelErr    = np.mean(RelErr_); # Moyenne des erreurs relatives
    titre_err    = "%s %s, ErG (Michel) : mean=%f min=%f, max=%f, std=%f" \
                    %(strset, varlib, MoyRelErr, np.min(RelErr_), \
                    np.max(RelErr_), np.std(RelErr_));
    print(titre_err);
    
    # Les valeurs par image (Daily sum)
    Somimg_ref = np.sum(phi_ref, axis=(1,2));
    Somimg_est = np.sum(phi_est, axis=(1,2));
    MatCor     = np.corrcoef(Somimg_ref, Somimg_est);
    correl     = MatCor[0,1];
    titre_valimg = "%s %s: Daily sum ; Coef. Correl. Ref.vs Est. : %f"%(strset, varlib, correl);
    print(titre_valimg);

    if FIGBYIMGS : #! Par image ordre chrono        
        # Par image dans l'ordre chronologique
        ids = np.argsort(calX[1]); # indice de tri des dates 
        plt.figure();
        plt.plot(np.array(RelErr_)[ids], '-*');
        plt.plot([0,Nimg_-1],[MoyRelErr, MoyRelErr]); # Trait de la moyenne
        plt.xticks(np.arange(Nimg_),calX[0][ids], rotation=35, ha='right');
        plt.axis("tight");
        plt.title(titre_err);
        #
        # Les valeurs par image
        plt.figure();        
        plt.plot(Somimg_ref[ids], '-*');        
        plt.plot(Somimg_est[ids], '-*');
        plt.xticks(np.arange(Nimg_),calX[0][ids], rotation=35, ha='right');
        plt.axis("tight");
        plt.legend(["NATL60", "Resac"], framealpha=0.5);
        plt.title(titre_valimg);
  
    # Somme Temporelle pour chaque pixel (pour une �tude spatiale)
    if not FIGBYPASS and reso_ in ResoAll : 
        PHI_ = np.sum(phi_ref, axis=(0));
        plt.figure(); plotavar(PHI_, "", None, None, None, cmap=wcmap[wk_]);
        plt.title("%s NATL60: Somme temporelle pour l'%s des Pixels"%(strset,varlib));            
        PHI_ = np.sum(phi_est, axis=(0));
        plt.figure(); plotavar(PHI_, "", None, None, None, cmap=wcmap[wk_]);
        plt.title("%s Resac: Somme temporelle pour l'%s des Pixels"%(strset,varlib));

    # Pourcentatge au Nord et au Sud d'une Latitude donn�e
    if LIMLAT_NORDSUD > 0 :             
        NL_lim   = nl2limlat(Nlig_, LIMLAT_NORDSUD); # Nombre de ligne jusqu'� la latitude limite Nord-Sud
        # for NATL60
        sumphi_ref = np.sum(phi_ref);
        XN_, XS_   = splitns(phi_ref, NL_lim);
        sumS       = np.sum(XS_);           sumN   = np.sum(XN_);
        pcentS     = sumS / sumphi_ref;     pcentN = sumN / sumphi_ref;
        print("%s Pcent %s NATL60: Nord=%.4f ; Sud=%.4f ; (sum pcent = %.4f (should be 1.)"
              %(strset, varlib, pcentN, pcentS, pcentN+pcentS));
        # for RESAC
        sumphi_est = np.sum(phi_est);
        XN_, XS_   = splitns(phi_est, NL_lim);            
        sumS       = np.sum(XS_);           sumN   = np.sum(XN_);
        pcentS     = sumS / sumphi_est;     pcentN = sumN / sumphi_est;
        print("%s Pcent %s RESAC : Nord=%.4f ; Sud=%.4f ; (sum pcent = %.4f (should be 1.)"
              %(strset, varlib, pcentN, pcentS, pcentN+pcentS));
        del XN_, XS_;               
   
    # Pour Histogramme "compar�s" (Ref., Est.) et scatter plot
    if flag_histo or flag_scat :
        Yref_ = phi_ref.ravel()
        Yest_ = phi_est.ravel()
        
    if flag_histo : # Histogramme "compar�s" (Ref., Est.)
        epsi = 0.0;
        min_ = np.min([phi_ref, phi_est]);
        print("Histo %s Ref-Est : min = %f"%(varlib, min_));
        if min_ <= 0 :
            epsi = -min_ + 1e-30;
        # histo en log ...                                                             
        logrefe = np.log(Yref_ + epsi);
        logeste = np.log(Yest_ + epsi);
        histodbl(logrefe, logeste, legende=['Ref','Est'], nbins=NBINS);
        plt.title("%s Histo log(%s+epsi=%e) (%dbins),\nrelerr_mean=%.4e"
                  %(strset,varlib, epsi, nbins, MoyRelErr), fontsize=x_figtitlesize);
    #
    if flag_scat : # Scatter plot not in log 
        pipo, pipo, R2 = scatplot(Yref_, Yest_)
        plt.title("%s %s scatter-plot, relerr_mean=%.4e R2=%f"%(strset,varlib,MoyRelErr,R2), fontsize=x_figtitlesize);               
#----------------------------------------------------------------------
def anycorrplex (Urefout, Vrefout, Uestout, Vestout, strset, calX=None) :
    # Elargissement de la Corr�lation Complexe � d'autre variable que u et v ...
    Nimg_ = len(Urefout);
    def corrplexpo (U, Up) :
        import cmath
        U_  = U.ravel();
        Up_ = Up.ravel()
        Upc = Up_.conjugate(); # Up conjug�
        Cc  = np.dot(U_,Upc) / (np.linalg.norm(U_) * np.linalg.norm(Up_)); # Corr�lation complexe
        modulCc, argCc = cmath.polar(Cc); # module->corr�lation, argCc->angle  \ entre les 2 champs complexes
        return modulCc, argCc, Cc, Upc;
    U_  = Urefout + 1j*Vrefout;
    Up_ = Uestout + 1j*Vestout;
    modulCc, argCc, Cc, Upc = corrplexpo (U_, Up_);
            
    # Calcul par image, always
    ModCC = []; ACC = [];
    for k in np.arange(Nimg_) : # Pour chaque image k           
        Uk_  = Urefout[k] + 1j*Vrefout[k];
        Upk_ = Uestout[k] + 1j*Vestout[k];
        modulCck, argCck, Cck, Upck = corrplexpo (Uk_, Upk_);
        ModCC.append(modulCck);
        ACC.append(argCck);
    # Modules
    moymodcc = np.mean(ModCC);
    titre_mod = "%s Correlation complexe: Module : Moy=%f, Min=%f, Max=%f, Std=%f" \
                %(strset, moymodcc, np.min(ModCC), np.max(ModCC), np.std(ModCC));
    print(titre_mod);
    # Angles en radian (no more)          
    # Angles en degr� absolu
    ACC       = np.array(ACC);
    AbsACC    = np.abs(ACC*180/np.pi); 
    moyabsacc = np.mean(AbsACC);
    titre_ang = "%s Correlation complexe: |Angle(degre)| : Moy=%f, Min=%f, Max=%f, Std=%f" \
                %(strset, moyabsacc, np.min(AbsACC),np.max(AbsACC),np.std(AbsACC));
    print(titre_ang);

    if FIGBYIMGS : # Par image
        # Par image ordre chrono
        ids = np.argsort(calX[1]); # indice de tri des dates
        # Module
        ModCC = np.array(ModCC)[ids];
        plt.figure(); plt.plot(ModCC,'.-');          
        plt.plot([0,Nimg_-1],[moymodcc, moymodcc]); 
        plt.axis("tight");
        plt.xticks(np.arange(Nimg_),calX[0][ids], rotation=35, ha='right');       
        plt.title(titre_mod);
        
        # Angles en degr� absolu 
        AbsACC = AbsACC[ids]
        plt.figure(); plt.plot(AbsACC,'.-');
        plt.plot([0,Nimg_-1],[moyabsacc, moyabsacc]); 
        plt.axis("tight");
        plt.xticks(np.arange(Nimg_),calX[0][ids], rotation=35, ha='right');
        plt.title(titre_ang);
#----------------------------------------------------------------------
def UVresult(Urefout, Vrefout, Uestout, Vestout, strset, flag_cosim=0, 
             flag_nrj=0, flag_ens=0, flag_corrplex=0, flag_div=0, dx=dxRout, dy=dyRout,
             flag_scat=0, flag_histo=0, calX=None, wdif='log', nbins=NBINS) :
    print("-- in UVresult --"); #print(np.shape(Urefout)); #(55L, 1L, 144L, 153L)
    if len(np.shape(Urefout)) == 4 :
        Nimg_, Ncan_, Nlig_, Ncol_ = np.shape(Urefout); 
    else : #== 3 :
        Nimg_, Nlig_, Ncol_ = np.shape(Urefout);
    reso_ = NLigR01/Nlig_;
    Npix_ = Nlig_*Ncol_;
    Nall_ = Nimg_*Npix_; 
    
    if flag_cosim :# Cosine Similarity
        dcos = 1 - dcosine2D(Urefout, Vrefout, Uestout, Vestout);
        print("%s cosine similarity on brute scalled : %.6f"%(strset,dcos));

    if flag_nrj : # Energie
        nrj_ref = (Urefout**2 + Vrefout**2) / 2; print("nrj shape : ", np.shape(nrj_ref));
        nrj_est = (Uestout**2 + Vestout**2) / 2;
        phistuff (nrj_ref, nrj_est, "NRJ", "Energie", strset,
                  calX, wdif, flag_histo, flag_scat, nbins=nbins);

    if flag_ens : # Enstrophie
        enst_ref = enstrophie2d(Urefout, Vrefout, dx, dy); #print("enstro shape : ", np.shape(enst_ref));
        enst_est = enstrophie2d(Uestout, Vestout, dx, dy);
        phistuff(enst_ref, enst_est, "ENS", "Enstrophie", strset,
                  calX, wdif, flag_histo, flag_scat, nbins=nbins);

    if flag_nrj and flag_ens and FIGBYIMGS :
        # Correlation Energie \ Enstrophie
        SomEimg_ref = np.sum(nrj_ref, axis=(1,2,3));
        SomEimg_est = np.sum(nrj_est, axis=(1,2,3));
        SomSimg_ref = np.sum(enst_ref, axis=(1,2));
        SomSimg_est = np.sum(enst_est, axis=(1,2));
        r_ref_ = np.corrcoef(SomEimg_ref, SomSimg_ref); 
        r_est_ = np.corrcoef(SomEimg_est, SomSimg_est);
        print("Correlation temporelle Energie \ Enstrophie : NATL60:=%.2f ; Resac:=%.2f"
              %(r_ref_[0,1],r_est_[0,1]))

    if flag_div : # Divergence
        div_ref    = divergence2d(Urefout, Vrefout, dx, dy);        
        div_est    = divergence2d(Uestout, Vestout, dx, dy);
        RMS, Nnan, inan = nanrms(div_ref, div_est);
        print("%s Divergence RMS (%dnan) : %.4e" %(strset, Nnan, RMS));

    if flag_corrplex : # Correlation Complexe
        anycorrplex (Urefout, Vrefout, Uestout, Vestout, strset, calX=calX); 
#
#--------------------------------------------------
def resultuv(Yrefout, Yestout, strset, flag_cosim=0, flag_nrj=0,  
             flag_ens=0, flag_corrplex=0, flag_div=0, dx=dxRout, dy=dyRout,
             flag_scat=0, flag_histo=0, calX=None, wdif='log', nbins=NBINS) :
    print("-- in resultuv --");
    iu_ = varOut.index("SSU");  Urefout = Yrefout[iu_];  Uestout = Yestout[iu_];
    iv_ = varOut.index("SSV");  Vrefout = Yrefout[iv_];  Vestout = Yestout[iv_];       
    if LIMLAT_NORDSUD <= 0 :
        UVresult(Urefout, Vrefout, Uestout, Vestout, strset, flag_cosim=flag_cosim,
                 flag_nrj=flag_nrj, flag_ens=flag_ens,
                 flag_corrplex=flag_corrplex, flag_div=flag_div, dx=dx, dy=dy,
                 flag_scat=flag_scat, flag_histo=flag_histo, calX=calX,
                 wdif=wdif, nbins=nbins);
    else : # PLM, je choisis de faire que tout, ou que Nord-Sud, sinon on
           # croule sous les images
        pipo, pipo, Nlig_, pipo = np.shape(Urefout);
        # Nombre de ligne jusqu'� la latitude limite Nord-Sud
        NL_lim = nl2limlat (Nlig_, LIMLAT_NORDSUD);
        # Split Nord-Sud
        UrefN_,  UrefS_ = splitns(Urefout, NL_lim);
        VrefN_,  VrefS_ = splitns(Vrefout, NL_lim);
        UestN_,  UestS_ = splitns(Uestout, NL_lim);
        VestN_,  VestS_ = splitns(Vestout, NL_lim);
        # Then
        #strset = strset+"-Nord";
        UVresult(UrefN_, VrefN_, UestN_, VestN_, strset+"-Nord", flag_cosim=flag_cosim,
             flag_nrj=flag_nrj, flag_ens=flag_ens,
             flag_corrplex=flag_corrplex, flag_div=flag_div, dx=dx, dy=dy,
             flag_scat=flag_scat, flag_histo=flag_histo, calX=calX,
             wdif=wdif, nbins=nbins);
        #strset = strset+"-Sud";
        UVresult(UrefS_, VrefS_, UestS_, VestS_, strset+"-Sud", flag_cosim=flag_cosim,
             flag_nrj=flag_nrj, flag_ens=flag_ens,
             flag_corrplex=flag_corrplex, flag_div=flag_div, dx=dx, dy=dy,
             flag_scat=flag_scat, flag_histo=flag_histo, calX=calX,
             wdif=wdif, nbins=nbins); 
#**********************************************************************
#----------------------------------------------------------------------
#
