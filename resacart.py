# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import sys
from   time  import time
import datetime
import numpy as np
import matplotlib.pyplot as plt
from   matplotlib  import cm
from   resacartparm  import *
from   resacartdef import *
#======================================================================
#######################################################################
#           GET AND SET THE REQUIRED BRUTE DATA
#######################################################################
#======================================================================
#%% Lecture Des Donn�es
Data_       = np.load("../donnees/natl60_htuv_01102012_01102013.npz") # A changer en fonction du chemin d'accès aux données
FdataAllVar = Data_['FdataAllVar']
varlue      = list(Data_['varlue'])
varlue      = [i.decode() for i in varlue] # pour éviter qu'il y ait un b devant les str
varlue[2]="U"
varlue[3]="V"

Nvar_, Nimg_, Nlig_, Ncol_ = np.shape(FdataAllVar) #(4L, 366L, 1296L, 1377L) 
#%%          


if FLAG_STAT_BASE_BRUTE : # Statistiques de base sur les variables brutes lues (->PPT)
    print("Statistique de base sur les variables brutes lues")
    stat2base(FdataAllVar, varlue)
    if LIMLAT_NORDSUD > 0 : # stat de base Nord-Sud
        for i in np.arange(len(FdataAllVar)) :
            XiNord_, XiSud_ = splitns(FdataAllVar[i], LIMLAT_NORDSUD)
            print("%s Nord: "%varlue[i], end='') 
            statibase(XiNord_)
            print("%s Sud : "%varlue[i], end='') 
            statibase(XiSud_)
        del XiNord_, XiSud_
#%%
# Splitset Ens App - Val - Test
print("Splitset Ens App - Val - Test ...", end='')
indA, indV, indT = isetalea(Nimg_, pcentSet)
VA_brute = [] 
VV_brute = [] 
VT_brute = []
for i in np.arange(Nvar_) : # Pour chaque variable (i.e. liste) (dans l'ordre de tvwmm ...)
    VA_brute.append(FdataAllVar[i][indA])              
    VV_brute.append(FdataAllVar[i][indV])              
    VT_brute.append(FdataAllVar[i][indT])
calA_ = calendrier[indA]
calA = [calA_, indA]
calV_ = calendrier[indV]
calV = [calV_, indV]
calT_ = calendrier[indT]
calT = [calT_, indT]
print("done (and for calendar too)")
#
del FdataAllVar #<<<<<<<
#
if FLAG_STAT_BASE_BRUTE_SET :
    # Stats de base en donn�es brute par ensemble
    print("APP :") 
    stat2base(VA_brute, varlue)
    print("VAL :")
    stat2base(VV_brute, varlue)
    print("TEST:")
    stat2base(VT_brute, varlue)
#
if FLAG_HISTO_VAR_BRUTE_SET :
    # Hist de comparaison des distributions par variable et par ensemble (->draft, ppt)
    for i in np.arange(len(varlue)) :
        plt.figure() 
        plt.suptitle("Histo %s"%varlue[i])
        
        H_ = VA_brute[i]
        plt.subplot(3,1,1) 
        plt.hist(H_.ravel(), bins=50)
        plt.title("APP")
        
        H_ = VV_brute[i]
        plt.subplot(3,1,2)
        plt.hist(H_.ravel(), bins=50)
        plt.title("VAL")
        
        H_ = VT_brute[i]
        plt.subplot(3,1,3) 
        plt.hist(H_.ravel(), bins=50)
        plt.title("TEST")
    plt.show()
#
# Make resolution for IN and OUT
VAout_brute, VVout_brute, VTout_brute, VAin_brute, VVin_brute, VTin_brute \
= setresolution(VA_brute,VV_brute,VT_brute,varlue,ResoIn,ResoOut)

# NATL60 : Moyennes Journali�res des pixels pour l'energie et
# l'enstrophie � la r�solution de sortie (R09)
if FLAG_DAILY_MOY_EE and IS_UVout : 
    iu_ = varOut.index("U")
    iv_ = varOut.index("V")
    Uall_ = np.concatenate((VAout_brute[iu_],VVout_brute[iu_],VTout_brute[iu_])) 
    Vall_ = np.concatenate((VAout_brute[iv_],VVout_brute[iv_],VTout_brute[iv_])) 
    NI_, NL_, NC_ = np.shape(Uall_)
    NL_lim_ = nl2limlat(NL_, LIMLAT_NORDSUD) # Nombre de ligne jusqu'� la latitude limite Nord-Sud
    #1) Energie
    NRJ_    = (Uall_**2 + Vall_**2) / 2
    moyNRJ_ = np.mean(NRJ_, axis=0)
    plt.figure() 
    plotavar(moyNRJ_, "", None, None, None)
    plt.title("NATL60: Daily Mean of Pixels Energy (at Rout=R09)")
    # Pourcentatge au Nord et au Sud d'une Latitude donn�e
    XN_, XS_ = splitns(NRJ_, NL_lim_)
    sumNRJ_  = np.sum(NRJ_)
    sumNRJS_ = np.sum(XS_)
    sumNRJN_ = np.sum(XN_)
    pcentS_  = sumNRJS_ / sumNRJ_     
    pcentN_  = sumNRJN_ /sumNRJ_
    print("Pcent Energy NATL60 (at Rout=R09): Nord=%.4f  Sud=%.4f  (sum pcent = %.4f (should be 1.)"
          %(pcentN_, pcentS_, pcentN_+pcentS_))
    del NRJ_, moyNRJ_, XN_, XS_, sumNRJ_, sumNRJS_, sumNRJN_   
    #
    #2) Enstrophie
    ENS_    = enstrophie2d(Uall_, Vall_, dxRout, dyRout)
    moyENS_ = np.mean(ENS_, axis=0)
    plt.figure() 
    plotavar(moyENS_, "", None, None, None)
    plt.title("NATL60: Daily Mean of Pixels Enstrophy (at Rout=R09)")
    # Pourcentatge au Nord et au Sud d'une Latitude donn�e
    XN_, XS_ = splitns(ENS_, NL_lim_)
    sumENS_  = np.sum(ENS_)
    sumENSS_ = np.sum(XS_)
    sumENSN_ = np.sum(XN_)
    pcentS_  = sumENSS_ / sumENS_
    pcentN_  = sumENSN_ /sumENS_
    print("Pcent Enstrophy NATL60 (at Rout=R09) : Nord=%.4f  Sud=%.4f  (sum pcent = %.4f (should be 1.)"
          %(pcentN_, pcentS_, pcentN_+pcentS_))
    del ENS_, moyENS_, XN_, XS_, sumENS_, sumENSS_, sumENSN_
    del Uall_, Vall_
#%%
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
print("# Mise en forme") 
for i in np.arange(NvarIn) :  # App In
    NdonA, pixlinA, pixcinA  = np.shape(VAin_brute[i])   # Nombre de donn�es et tailles
    VAin_brute[i] = VAin_brute[i].reshape(NdonA,1,pixlinA,pixcinA)
for i in np.arange(NvarOut) : # App Out
    NoutA, pixloutA, pixcoutA = np.shape(VAout_brute[i])
    VAout_brute[i] = VAout_brute[i].reshape(NdonA,1,pixloutA,pixcoutA)   
if NdonA != NoutA :
    raise ValueError("Probl�me A") # ce n'est pas suffisant

if TEST_ON :
    for i in np.arange(NvarIn) :  # Tst In
        NdonT, pixlinT, pixcinT  = np.shape(VTin_brute[i])     # Nombre de donn�es et tailles
        VTin_brute[i] = VTin_brute[i].reshape(NdonT,1,pixlinT,pixcinT)
    for i in np.arange(NvarOut) : # Tst Out       
        NoutT, pixloutT, pixcoutT = np.shape(VTout_brute[i])
        VTout_brute[i] = VTout_brute[i].reshape(NdonT,1,pixloutT,pixcoutT)
    if NdonT != NoutT :
        raise ValueError("Probl�me T") # ce n'est pas suffisant
if VALID_ON : 
    for i in np.arange(NvarIn) :  # Val In
        NdonV, pixlinV, pixcinV  = np.shape(VVin_brute[i]) # Nombre de donn�es et tailles
        VVin_brute[i] = VVin_brute[i].reshape(NdonV,1,pixlinV,pixcinV)
    for i in np.arange(NvarOut) : # Val Out       
        NoutV, pixloutV, pixcoutV = np.shape(VVout_brute[i])
        VVout_brute[i] = VVout_brute[i].reshape(NdonV,1,pixloutV,pixcoutV)
    if NdonV != NoutV :
        raise ValueError("Probl�me V") # ce n'est pas suffisant
#
#======================================================================
# Ajout de bruit sur l'input SSH brute pour l'ensemble de TEST uniquement
if SIGT_NOISE > 0 :
    for i in np.arange(NvarIn) :
        if varIn[i] == "SSH" :
            STAT_ON_NOISE = True
            if STAT_ON_NOISE : 
                X0_ = VTin_brute[i]  # H avant bruitage
            #
            np.random.seed(0)
            N_,M_,L_,C_ = np.shape(VTin_brute[i])
            #sigma = 0.050 #1.=8�, 10=42�, 100=45�(0�-180�)
            noise_ = np.random.randn(N_,M_,L_,C_) * SIGT_NOISE
            VTin_brute[i] = VTin_brute[i] + noise_     # TEST set only
            #
            if STAT_ON_NOISE :            
                noiseabs_ = np.abs(noise_)
                print("SSH : noise %.3f (min=%.3f, max=%.3f, moy=%.4e, std=%.3f moyabs=%.4e, stdabs=%.3f)"
                      %(SIGT_NOISE, np.min(noise_), np.max(noise_), np.mean(noise_),
                        np.std(noise_), np.mean(noiseabs_), np.std(noiseabs_)))
                rmsi, Nnan, inan = nanrms(VTin_brute[i], X0_)
                print("noiseRMSE : %f "%(rmsi))
                #          
                print("noiseErrRelAbs : %f "%( np.sum(noiseabs_) / np.sum(np.abs(X0_))))
                print("Stat noiseAbs perturbation (sur H en m): min=%.4e, max=%.4f, moy=%.4f, std=%.4f"
                      %(np.min(noiseabs_), np.max(noiseabs_),np.mean(noiseabs_), np.std(noiseabs_)))
                del noiseabs_, X0_
            del noise_ 
#======================================================================
if VisuBA+VisuBV+VisuBT > 0 :
    print("# Visualisation des donn�es brutes")
    
#%%
#----------------------------------------------------------------------
def visuB (X_brute, varIO, VisuB, Ndon, strset, inout, im2show,
           qmask=None, qscale=None, qmode=None, calX0=None) :
    Nvar = len(varIO)
    calX0_= None
    if calX0 is not None :
        calX0_ = calX0[im2show]
    if (VisuB==1 or VisuB==3) :
        for i in np.arange(Nvar) :                
            wk_ = tvwmm.index(varIO[i])
            showimgdata(X_brute[i],cmap=wcmap[wk_], n=Ndon,fr=0, vmin=wbmin[wk_],
                        vmax=wbmax[wk_], vnorm=wnorm[wk_], origine='lower')            
            plt.suptitle("%s: %s_brute(%s%d), min=%f, max=%f, mean=%f, std=%f"
                    %(strset, varIO[i], inout, i+1, np.min(X_brute[i]),np.max(X_brute[i]),
                      np.mean(X_brute[i]),np.std(X_brute[i])),
                      fontsize=x_figtitlesize)
   
    if (VisuB==2 or VisuB==3) and len(im2show) > 0 :
        if inout=="in" or (inout=="out" and 1) :
            for i in np.arange(Nvar) :             
                wk_ = tvwmm.index(varIO[i])
                suptitre="some %s brute(%s%d) %s"%(strset,inout,i+1,im2show)
                showsome(X_brute[i][im2show,0,:,:], wmin=wbmin[wk_], wmax=wbmax[wk_],
                         wnorm=wnorm[wk_], cmap=wcmap[wk_], fsize=(12,7), calX0=calX0_,
                         varlib=varIO[i], suptitre=suptitre, Xtit="")
                
        if (VisuB==2 or VisuB==3) and len(im2show) > 0  and inout=="out" :     
            if IS_HUVout and 1 : # Vecteurs UV on H       
                ih_ = showhquiv(X_brute, im2show, qscale=qscale, qmask=qmask, qmode=qmode, calX0=calX0)    
                plt.suptitle("some %s(out) %s_brute(+uv) %s"%(strset, varOut[ih_], im2show),
                         fontsize=x_figtitlesize)   
#--------
if VisuBA > 0 :
    # visu des donn�es brutes d'APP en ENTREE
    visuB(VAin_brute, varIn, VisuBA, NdonA, "APP","in", im2showA, calX0=calA[0])
    # Visu des donn�es brutes d'APP en SORTIE.
    visuB(VAout_brute, varOut, VisuBA, NdonA, "APP","out", im2showA, calX0=calA[0]) 
if VisuBT > 0 and TEST_ON : 
    # visu des donn�es brutes de TEST en entr�e
    visuB(VTin_brute, varIn, VisuBT, NdonT, "TEST","in", im2showT, calX0=calT[0]) 
    # Visu des donn�es brutes de TEST en sortie.
    visuB(VTout_brute, varOut, VisuBT, NdonT, "TEST","out", im2showT, calX0=calT[0]) 
if VisuBV > 0 and VALID_ON : 
    # visu des donn�es brutes de VALIDATION en entr�e
    visuB(VVin_brute, varIn, VisuBV, NdonV, "VALID","in", im2showV, calX0=calV[0]) 
    # Visu des donn�es brutes de VALIDATION en sortie.
    visuB(VVout_brute, varOut, VisuBV, NdonV, "VALID","out", im2showV, calX0=calV[0]) 
if VisuBstop or SCENARCHI==0 :
    print("STOP after visu donn�es brutes")
    plt.show() 
    sys.exit(0)
#%%
#======================================================================
#                   CODIFICATION / NORMALISATION
#======================================================================
print("# Codification / Normalisation")
# PLM, CE DOIT ETRE OBLIGATOIRE car la sauvegarde des param�tres n'est 
# pas faite, Il faut repasser ici pour les recalculer � chaque fois
VAin = []
coparmAin = []
for i in np.arange(NvarIn) :
    VAin_,  coparmAin_  = codage(VAin_brute[i],  "fit01")
    print(coparmAin_)
    VAin.append(VAin_)
    coparmAin.append(coparmAin_)
del VAin_, coparmAin_
x_train = VAin
NcanIn  = len(x_train)
#
VAout = [] 
coparmAout = []
for i in np.arange(NvarOut) :
    VAout_,  coparmAout_  = codage(VAout_brute[i],  "fit01") 
    print(coparmAout_)
    VAout.append(VAout_)
    coparmAout.append(coparmAout_)
del VAout_, coparmAout_  
y_train = VAout
NensA   = len(y_train[0])
#
if TEST_ON : # Il faut appliquer le m�me codage et dans les m�mes conditions
    # (i.e. avec les m�mes param�tres) que ceux de l'apprentissage.                                                                                 
    VTin = []
    for i in np.arange(NvarIn) :
        VTin_ = recodage(VTin_brute[i], coparmAin[i])
        VTin.append(VTin_)
    del VTin_
    x_test = VTin
    #
    VTout = []
    for i in np.arange(NvarOut) :
        VTout_ =  recodage(VTout_brute[i], coparmAout[i])
        VTout.append(VTout_)
    del VTout_
    #    
    y_test = VTout
    NensT = len(y_test[0])
    
if VALID_ON : # Il faut appliquer le m�me codage et dans les m�mes conditions
    # (i.e. avec les m�mes param�tre) que ceux de l'apprntissage.                                                                                  
    VVin = []
    for i in np.arange(NvarIn) :
        VVin_ = recodage(VVin_brute[i], coparmAin[i])
        VVin.append(VVin_)
    del VVin_
    x_valid = VVin
    #
    VVout = []
    for i in np.arange(NvarOut) :
        VVout_ =  recodage(VVout_brute[i], coparmAout[i])
        VVout.append(VVout_)
    del VVout_
    #
    y_valid = VVout
    NensV   = len(y_valid[0])
#----------------------------------------------------------------------
# POUR AVOIR CHANEL LAST, en Linux dans ~/.keras/keras.json
# Windowd c:/Users/charles/.keras/keras.json
#y_test_brute=[]
#y_train_brute=[]
for i in np.arange(NvarIn):
    x_train[i] = x_train[i].transpose(0,2,3,1)
    x_valid[i] = x_valid[i].transpose(0,2,3,1)
    x_test[i] = x_test[i].transpose(0,2,3,1)
for i in np.arange(NvarOut):
    y_train[i] = y_train[i].transpose(0,2,3,1)
    y_valid[i] = y_valid[i].transpose(0,2,3,1)
    y_test[i] = y_test[i].transpose(0,2,3,1)
#    y_test_brute.append( VTout_brute[i].transpose(0,2,3,1))
#    y_train_brute.append(VAout_brute[i].transpose(0,2,3,1))
#----------------------------------------------------------------------

if 1 : # Affichage des shapes ...
    for i in np.arange(NvarIn) : 
        print("%s shape x_train : "%varIn[i], np.shape(x_train[i]))
    for i in np.arange(NvarOut) : 
        print("%s shape y_train : "%varOut[i], np.shape(y_train[i]))
    if VALID_ON :
        for i in np.arange(NvarIn) : 
            print("%s shape x_valid : "%varIn[i], np.shape(x_valid[i]))
        for i in np.arange(NvarOut) : 
            print("%s shape y_valid : "%varOut[i], np.shape(y_valid[i]))
    if TEST_ON :
        for i in np.arange(NvarIn) : 
            print("%s shape x_test : "%varIn[i], np.shape(x_test[i]))
        for i in np.arange(NvarOut) : 
            print("%s shape y_test : "%varOut[i], np.shape(y_test[i]))
            
            
#%%
#======================================================================
#######################################################################
#                       THE ARCHITECTURE
#######################################################################
#======================================================================
print("# Build and compile Architecture")
from keras.layers import Input, Dense, Flatten, Reshape, AveragePooling2D, Dropout 
from keras.layers import Conv2D, MaxPooling2D, UpSampling2D #, Deconvolution2D
from keras.layers import Concatenate, concatenate 
from keras.models    import Model
from keras.models    import load_model
import keras.callbacks
from keras.callbacks import EarlyStopping, ModelCheckpoint
#----------------------------------------------------------------------
if 1 : # Avoir sa propre fonction d'activation
    from keras.layers import Activation
    from keras import backend as K
    from keras.utils.generic_utils import get_custom_objects
    def sig010(x) : # sigmo�d dans l'intervalle [1.10 , 0.10]
        return  (K.sigmoid(x) * 1.20) - 0.10
    get_custom_objects().update({'sig01': Activation(sig010)})
#----------------------------------------------------------------------
# D�termination de la prise en compte (ou pas) de la SST en entr�e
# de l'archi selon les r�solutions indiqu�es (dans resacparm.py).
all_Kinput_img = []
IS_SSTR81 = IS_SSTR27 = IS_SSTR09 = IS_SSTR03 =IS_SSTR01= False
for ii in np.arange(NvarIn) :
    #Ncan_, NL_, NC_ = np.shape(x_train[ii][0]) # CHANEL FIRST
    #all_Kinput_img.append(Input(shape=(Ncan_,NL_,NC_)))
    NL_, NC_, Ncan_  = np.shape(x_train[ii][0]) # CHANEL LAST
    all_Kinput_img.append(Input(shape=(NL_,NC_,Ncan_)))
    
    if varIn[ii]=="SST" :
        if ResoIn[ii]==81 :
            if IS_SSTR81 == False :
                ISSTR81   = all_Kinput_img[ii]
                IS_SSTR81 = True
            else :
                ISSTR81 = concatenate([ISSTR81, all_Kinput_img[ii]], axis=1)                   
        elif ResoIn[ii]==27 :
            if IS_SSTR27 == False :
                ISSTR27   = all_Kinput_img[ii]
                IS_SSTR27 = True
            else :
                ISSTR27 = concatenate([ISSTR27, all_Kinput_img[ii]], axis=1)
        elif ResoIn[ii]==9 :
            if IS_SSTR09 == False :
                ISSTR09   = all_Kinput_img[ii]
                IS_SSTR09 = True
            else :
                ISSTR09 = concatenate([ISSTR09, all_Kinput_img[ii]], axis=1)
        elif ResoIn[ii]==3 :
            if IS_SSTR03 == False :
                ISSTR03   = all_Kinput_img[ii]
                IS_SSTR03 = True
            else :
                ISSTR03 = concatenate([ISSTR03, all_Kinput_img[ii]], axis=1)
        #elif ResoIn[ii]==1:
        #    if IS_STR01 == False :
        #        ISSTR01 = all_Kinput_img[ii]
        #        IS_SSTR01 = True
        #    else:
        #        ISSTR01= concatenate([ISSTR01, all_Kinput_img[ii]], axis=1)
        else :
            raise ValueError("Other resolution of SST not prevue")
#
# Param�tres par d�faut (a re-adapter si besoin)
np.random.seed(acide)      # Choose a random (or not) for reproductibilitie
init       = 'he_normal'   #'glorot_normal'+ #'orthogonal' #'normal'
factiv     = 'relu'        # 'relu', 'tanh', 'sigmoid', 'linear'
factout    = sig010        # 'relu', 'tanh', 'sigmoid', 'linear'
upfactor   =  3            # facteur d'upsampling
#
ArchiOut   = []            # Liste des sorties de l'Archi
#
# 1er �tage : ... to SSH_R27
if IS_SSTR81 : # SSH_R81 + SST_R81 to SSH_R27
    ISSHreso = all_Kinput_img[0] # Input SSH reso
    #ArchiA   = concatenate([ISSHreso, ISSTR81], axis=1)# CHANEL FIRST
    ArchiA   = concatenate([ISSHreso, ISSTR81], axis=3)# CHANEL LAST
    ArchiA   = UpSampling2D((upfactor, upfactor))(ArchiA)
    for i in np.arange(9) : 
        ArchiA = Conv2D(36,(6,6), activation=factiv, padding='same',kernel_initializer=init)(ArchiA)
    ArchiA = Conv2D(1,(1,1), activation=factout, padding='same',kernel_initializer=init,name='archiA')(ArchiA)
elif IS_SSTR27: # # SSH_R81 + SST_R09 to SSH_R27
    ISSHreso = all_Kinput_img[0] # Input SSH reso  81
    ArchiA   = UpSampling2D((upfactor, upfactor))(ISSHreso)          
    #ArchiA   = concatenate([ArchiA, ISSTR27], axis=1);# CHANEL FIRST
    ArchiA   = concatenate([ArchiA, ISSTR27], axis=3)# CHANEL LAST
    for i in np.arange(7) : #ou (9) ?
        ArchiA = Conv2D(36,(6,6), activation=factiv, padding='same',kernel_initializer=init)(ArchiA)
    ArchiA = Conv2D(1,(2,2), activation=factout, padding='same',kernel_initializer=init,name='archi81to27')(ArchiA) 
else : # assumed Sans SST : # SSH_R81 to SSH_R27
    ISSHreso = all_Kinput_img[0] # Input SSH reso
    ArchiA   = UpSampling2D((upfactor, upfactor))(ISSHreso)
    for i in np.arange(7) :
        ArchiA = Conv2D(36,(6,6), activation=factiv, padding='same',kernel_initializer=init)(ArchiA)
    ArchiA = Conv2D(1,(1,1), activation=factout, padding='same',kernel_initializer=init,name='archiA')(ArchiA)
ArchiOut.append(ArchiA)

# 2�me �tage H_R27+SST_R09 to SSH_R09
ArchiB = UpSampling2D((upfactor, upfactor))(ArchiA)
if IS_SSTR09 and 1 : # SST R09 utilis� pour trouver SSH-R09
    #ArchiB   = concatenate([ArchiB, ISSTR09], axis=1);# CHANEL FIRST
    ArchiB   = concatenate([ArchiB, ISSTR09], axis=3)# CHANEL LAST
for i in np.arange(3) :
    ArchiB = Conv2D(24,(5, 5), activation=factiv, padding='same',kernel_initializer=init)(ArchiB)
ArchiB = Conv2D(1,(1, 1), activation=factout, padding='same',kernel_initializer=init,name='archi27to09')(ArchiB)
ArchiOut.append(ArchiB)


# 3eme etage SSH_R09+ SST_R03 to SSH_R03
ArchiC = UpSampling2D((upfactor,upfactor))(ArchiB)
if IS_SSTR03 and 1:
    ArchiC= concatenate([ArchiC,ISSTR03], axis=3)#CHANEL LAST
for i in np.arange(2):
    ArchiC= Conv2D(24,(5,5), activation=factiv, padding='same',kernel_initializer=init)(ArchiC)
ArchiC=Conv2D(12,(3,3), activation=factiv, padding='same',kernel_initializer=init)(ArchiC)
#((taille de l'image - taille kernel +2*padding)/ (stride)) + 1
ArchiC= Conv2D(1,(1,1), activation=factout, padding='same',kernel_initializer=init,name='archi09to03')(ArchiC)
ArchiOut.append(ArchiC)


# 4eme �tage : H_R03+SST_R03 to U_R03 and V_R03
nc, ml, mc, nl = (24, 5, 5, 2)
if IS_SSTR09 : # SST R09 utilis� pour trouver SSU R09 et SSV R09
    #ArchiUV  = concatenate([ArchiC, ISSTR03], axis=1);# CHANEL FIRST
    ArchiUV  = concatenate([ArchiC, ISSTR03], axis=3);# CHANEL LAST
    ArchiUV  = Conv2D(nc,(ml,mc), activation=factiv, padding='same',kernel_initializer=init)(ArchiUV) ##
else :
    ArchiUV = Conv2D(nc,(ml,mc), activation=factiv, padding='same',kernel_initializer=init)(ArchiC) ##
for i in np.arange(nl) :
    ArchiUV = Conv2D(nc,(ml,mc), activation=factiv, padding='same',kernel_initializer=init)(ArchiUV)
ArchiUV = Conv2D(1,(1, 1), activation=factout, padding='same',kernel_initializer=init)(ArchiUV)
#
ArchiU = Conv2D(nc,(ml,mc), activation=factiv, padding='same',kernel_initializer=init)(ArchiUV)
for i in np.arange(nl) :
    ArchiU = Conv2D(nc,(ml,mc), activation=factiv, padding='same',kernel_initializer=init)(ArchiU)
ArchiU = Conv2D(1,(1, 1), activation=factout, padding='same',kernel_initializer=init,name='archiU')(ArchiU)
ArchiOut.append(ArchiU)
#
ArchiV = Conv2D(nc,(ml,mc), activation=factiv, padding='same',kernel_initializer=init)(ArchiUV)
for i in np.arange(nl) :
    ArchiV = Conv2D(nc,(ml,mc), activation=factiv, padding='same',kernel_initializer=init)(ArchiV)
ArchiV = Conv2D(1,(1, 1), activation=factout, padding='same',kernel_initializer=init,name='archiV')(ArchiV); 
ArchiOut.append(ArchiV)
#
Mdl   = Model(all_Kinput_img, ArchiOut)
Mdl.summary(); 
Mdl.compile(loss='logcosh', optimizer='adam')
print("Architecture completed")
#%%
#
#======================================================================
#                LEARNING (ou reprendre)
print("\nLEARNING (ou reprendre) :", RUN_MODE)
#======================================================================
if RUN_MODE=="RESUME" :
    print("Reload des poids d'un model préalablement sauvegardé")
    Mdl.load_weights(Mdl2reload); 
    print("reload weights, done that")
#
elif RUN_MODE=="REPRENDRE":
    try:
        os.makedirs("REPRENDRE/Archi")
    except:
        print("Le dossier REPRENDRE/Archi existe deja")
    try:
        os.makedirs("REPRENDRE/Weights")
    except:
        print("Le dossier REPRENDRE/Weights existe deja")
    try:
        os.makedirs('REPRENDRE/History')
    except:
        print("Le dossier REPRENDRE/History existe deja")
    try:
        os.makedirs('REPRENDRE/Resultats')
    except:
        print("Le dossier REPRENDRE/Resultats existe deja")
    try:
        os.makedirs('REPRENDRE/Resultats/Images/')
    except:
        print("Le dossier REPRENDRE/Resultats/Images/ existe deja")
    print("Reprise d'un apprentissage passé")
    np.random.seed(acide)
    Mdl= load_model(Mdl2reprendre_archi)    # Chargement de l'archi
    Mdl.load_weights(Mdl2reprendre_poids)   # Chargement des poids à partir desquel on reprend l'apprentissage
    print("Chargement du modele a repreprendre effectue")
    with open('REPRENDRE/History/history.pkl', 'rb') as file: 
        HistoryFromReprendre=pickle.load(file)      #Chargement de l'historique des loss
    print("Chargement de l'historique des loss precedentes effectue")
    #
    
    t0=time()
    if VALID_ON==3: # Early stopping sur l'ensemle de validation. Les poids ou
                    # le mod�le sauvgardable sont obtenus � la fin du run
                    # apr�s patience it sans am�lioration).
        earlystop= EarlyStopping(monitor='val_loss',
                patience=int(Niter/4), verbose=2, mode='auto')
        H_reprendre = Mdl.fit(x_train, y_train, verbose=2, epochs=Niter,batch_size=Bsize,
                shuffle=True, callbacks=[earlystop], validation_data=(x_valid, y_valid))
    elif VALID_ON==4 or VALID_ON==5:
        log_dir = "REPRENDRE/Tensorboard/logs/fit/" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        tensorboard_callback = keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1)
        # 4 : Le run va jusqu'au bout. On r�cup�re, par la suite la sauvegarde
        #     des poids au meilleur de l'ensemble de validation pour les r�sultats.          
        checkpoint_reprendre= ModelCheckpoint("REPRENDRE/Resultats/"+"modelkvalid.ckpt", monitor='val_loss', verbose=2,
                save_best_only=True, save_weights_only=True, mode='auto')
        
        callback_history=LossHistory()
        try:
            os.remove(history_filename)
        except OSError:
            print("Le fichier n'existe pas encore il n'a pas pu être supprimé")
            pass
        
        callbacks_list=[checkpoint_reprendre, callback_history,tensorboard_callback] # Sauvegarde de H.history à chaque epoch
        
        H_reprendre = Mdl.fit(x_train, y_train, verbose=2, epochs=Niter, batch_size=Bsize,
                shuffle=True, callbacks=callbacks_list, validation_data=(x_valid, y_valid));          
    print("Le run du mode REPRENDRE est allé jusqu'au bout des %d itérations"%(Niter))
    print("learning time : %f secondes" %(time()-t0))
    print("Temps moyen par itération: %f secondes " %((time()-t0)/Niter))
    #
    print("La longueur du dictionnaire de l'historique de l'apprentissage 1 est ",len(HistoryFromReprendre['val_loss']))
    print("La longueur du dictionnaire de l'historique de l'apprentissage 2 est ",len(H_reprendre.history['val_loss']))
    temp_hist=HistoryFromReprendre
    for key in HistoryFromReprendre.keys():
        HistoryFromReprendre[key].extend(H_reprendre.history[key])
    
    print("#Visualisation de la courbes d'erreur en apprentissage et de performances")
    plt.figure()
    LH = []
    kh = HistoryFromReprendre.keys()
    LH.append(kh)
    if VALID_ON :
        valid_err_hist = HistoryFromReprendre['val_loss']
        itminval = np.argmin(valid_err_hist) 
        print("La val_loss n'a pas évolué après l'itération n° "+ str(itminval))
                
    for kle in HistoryFromReprendre.keys() :
        hloss = HistoryFromReprendre[kle]
        plt.plot(np.log(hloss),alpha=0.7)
        hloss = np.reshape(hloss,(1,len(hloss)))
        LH.append(hloss[0])               
    plt.legend(HistoryFromReprendre.keys(), framealpha=0.3, loc=0)
    plt.axvline(itminval,c="black")
    plt.xlabel('epoch')    
    plt.ylabel('log(loss)')
    plt.savefig("REPRENDRE/logloss_reprendre.png")
    if VALID_ON :
        valid_err_hist = HistoryFromReprendre['val_loss']
        itminval = np.argmin(valid_err_hist) 
        print("La val_loss n'a pas évolué après l'itération n° "+ str(itminval))
        
    #
    with open('REPRENDRE/Resultats/'+'history.pkl', 'wb') as file_pi:
        pickle.dump(HistoryFromReprendre, file_pi)
    if Mdl2reprendre_archi is not None:
        Mdl.save("REPRENDRE/Resultats/Archi_Model_From_Reprendre")#Sauvegarde du modèle après reprise
    
    


else : #=> RUN_MODE is "LEARN"
    print("#Learning");
    np.random.seed(acide);
    #
    try:
        os.makedirs(Mdl2save+"Historique_Loss/")
    except:
        print("Le dossier Historique_Loss existe deja")
    if Mdl2save is not None :
        # Sauvegarde de l'architecture du modele
        #Mdl.save_weights(Mdl2save+"wei");
        Mdl.save("Save_Model/Archi")
        print('Architecture du modèle sauvegardée')
        
        # Les param�tres de codage de l'ens d'App
        np.save(Mdl2save+"coparm", [coparmAin, coparmAout]);
        print("Paramètres du codage sauvegardés");
    t0 = time();    
    if VALID_ON == 0 : # Pas d'usage de l'ensemble de validation (run 'simple)
        H = Mdl.fit(x_train, y_train, verbose=2,
                    epochs=Niter,batch_size=Bsize,shuffle=True);
    elif VALID_ON == 1 : # Usage de l'ensemble de validation (sans early stopping)
        H = Mdl.fit(x_train, y_train, verbose=2, epochs=Niter, batch_size=Bsize,
                    shuffle=True, validation_data=(x_valid, y_valid));             
    elif VALID_ON == 3 : # Early stopping sur l'ensemle de validation. Les poids ou
                         # le mod�le sauvgardable sont obtenus � la fin du run
                         # apr�s patience it sans am�lioration).      
        earlystop = EarlyStopping(monitor='val_loss',
                patience=20, verbose=2, mode='auto');       
        H = Mdl.fit(x_train, y_train, verbose=2, epochs=Niter,batch_size=Bsize,
                shuffle=True, callbacks=[earlystop], validation_data=(x_valid, y_valid));
    elif VALID_ON==4 or VALID_ON==5 : 
        # 4 : Le run va jusqu'au bout. On r�cup�re, par la suite la sauvegarde
        #     des poids au meilleur de l'ensemble de validation pour les r�sultats.          
        log_dir = Mdl2save+"logs/fit/" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1)
        
        checkpoint= ModelCheckpoint(Mdl2save+"Weights/modelkvalid.ckpt", monitor='val_loss', verbose=2,
                save_best_only=True, save_weights_only=True, mode='auto')
        callback_history=LossHistory()
        
        try:
            os.remove(history_filename)
        except OSError:
            print("Le fichier n'existe pas encore il n'a pas pu être supprimé")
            pass
        
        callbacks_list=[checkpoint, callback_history] # Sauvegarde de H.history à chaque epoch
        H = Mdl.fit(x_train, y_train, verbose=2, epochs=Niter, batch_size=Bsize,
                shuffle=True, callbacks=callbacks_list, validation_data=(x_valid, y_valid));          
    print("learning time : %f" %(time()-t0))
    print("temps moyen des" +str(Niter)+ " iterations:"+str((time()-t0)/Niter))
    with open(Mdl2save+'history.pkl', 'wb') as file:
        pickle.dump(H.history, file)
    print("L'entrainement est allé jusqu'au bout des "+str(Niter)+" itérations" )
    #
    print("#Visualisation de la courbes d'erreur en apprentissage et de performances")
    plt.figure()
    LH = []
    kh = H.history.keys()
    LH.append(kh);            
    for kle in H.history.keys() :
        hloss = H.history[kle]
        plt.plot(np.log(hloss),alpha=0.8)
        hloss = np.reshape(hloss,(1,len(hloss)))
        LH.append(hloss[0])               
    plt.legend(H.history.keys(), framealpha=0.3, loc=0)
    ##np.save(Mdl2save+"loss",LH[1])#sauvegarde LOSS
    ##np.save(Mdl2save+"accuracy",LH[2])#sauvegarde accuracy
    #np.save(Mdl2save+"lhist", LH); # Sauvegarde list loss history
    plt.xlabel('epoch')    
    plt.ylabel('log(loss)')
    if VALID_ON :
        valid_err_hist = H.history['val_loss']
        itminval = np.argmin(valid_err_hist) 
        print("itminval =",itminval)
#%%
#----------------------------------------------------------------------
if RUN_MODE == "RESUME" :
    try : # loss history
        #LH = np.load(Mdl2reload+"lhist.npy", allow_pickle=True);
        with open('Save_Model/history.pkl', 'rb') as file:
            LH=pickle.load(file)
    except :
        print("failed when reloading file history to replot loss curves")
    try:
        os.makedirs("RESUME/Images")
    except:
        print("Le dossier RESUME/Images existe deja")
    plot_history_RESUME(LH)
    #plt.figure();
    #kh = list(LH.keys());
    #for i in np.arange(len(kh))+1 :
    #    hloss = LH[i];
    #    plt.plot(np.log(hloss),alpha=0.7);
    #plt.legend(kh, framealpha=0.3, loc=0);
    #plt.xlabel('epoch');    plt.ylabel('log(loss)');
    #ivloss = kh.index('val_loss');
    #vloss = LH[ivloss+1]
    #itminval = np.argmin(vloss); print("itminval =",itminval)
    
    try:
        os.makedirs("RESUME/Images")
    except:
        print("Le dossier RESUME/Images existe deja")

#%%

#======================================================================
#                       RESULTS ON MODEL LEARNED
#======================================================================
if VALID_ON == 4 :
    if RUN_MODE == "LEARN" :
        # ON prend les poids au min de la validation (sauvegard�s pendant le run)
        Mdl.load_weights(Mdl2save+"Weights/modelkvalid.ckpt");
    elif RUN_MODE == "REPRENDRE":
        #Si le resultat obtenu precedement est meilleur, on charge les poids precedents
        #Par contre si les poids après réentrainement sont meilleurs c'est eux que l'on charge
        if temp_hist['val_loss'][np.argmin(temp_hist['val_loss'])] < H_reprendre.history['val_loss'][np.argmin(H_reprendre.history['val_loss'])]:
            #Les anciens resultats sont meilleurs
            Mdl.load_weights(Mdl2reprendre_poids)
            print("On recharge les poids précédents car la reprise n'a pas amélioré la val_loss")
        else: # Le modele reprendre a eu de meilleures performances
            Mdl.load_weights("REPRENDRE/Resultats/modelkvalid.ckpt")
            print("On charge les nouveaux poids après reprise car la val_loss s'est améliorée")
    else : # suppose="RESUME" : on reprend des poids d'un run pr�cedent 
        #Mdl.load_weights(Mdl2reload+"valid.ckpt");
        print("Mode "+ RUN_MODE+" les poids ont été chargés")
    print("load weights min valid, done this");
#----------------------------------------------------------------------
# 
if RUN_MODE=="LEARN":  
    plot_history(H)

    print("Affichage des performances")
    #with open(Mdl2save+'test.pkl', 'wb') as file_pii:
    #    pickle.dump(x_test, file_pii)
    #print("Sauvegarde de la base de test réussie")
     
def setresult(Mdl, x_set, strset, VXout_brute, calX, flag_rms=0, flag_scat=0,
              flag_histo=0, nbins=NBINS, flag_cosim=0, flag_nrj=0, flag_ens=0,
              flag_corrplex=0, flag_div=0, dx=dxRout, dy=dyRout, wdif='log', im2show=[], savefig=False) :
    # Predicted coded data 
    y_scale = Mdl.predict(x_set);
    for i in np.arange(NvarOut):
        y_scale[i]=y_scale[i].transpose(0,3,1,2)
    if len(varOut)==1 : # Si une seule sortie mettre y_scale en forme de list comme l'est y_train
        y_scale = [y_scale];
    #BACK TO BRUTE

    print("%s BACK TO BRUTE"%strset);
 
    for i in targetout : # boucle sur les variables cibles en sortie            
        wk_ = tvwmm.index(varOut[i]);
        y_scale[i] = decodage(y_scale[i], coparmAout[i]);

        if not RESBYPASS :
            if LIMLAT_NORDSUD <= 0 :
                result("brute scalled", strset, varOut[i], VXout_brute[i], y_scale[i],
                       flag_rms=flag_rms, flag_scat=flag_scat, flag_histo=flag_histo,
                       nbins=nbins, calX=calX,savefig=savefig);
            else : # PLM, je choisis de faire que tout, ou que Nord-Sud, pour �viter
                   # de crouler sous les images
                pipo, pipo, Nlig_, pipo = np.shape(VXout_brute[i]);
                # Nombre de ligne jusqu'� la latitude limite Nord-Sud
                NL_lim = nl2limlat (Nlig_, LIMLAT_NORDSUD);
                # Split Nord-Sud
                youtN_,  youtS_ = splitns(VXout_brute[i], NL_lim); 
                yscalN_, yscalS_= splitns(y_scale[i], NL_lim);
                result("brute scalled", strset+"-Nord", varOut[i], youtN_, yscalN_,
                       flag_rms=flag_rms, flag_scat=flag_scat, flag_histo=flag_histo,
                       nbins=nbins, calX=calX);
                result("brute scalled", strset+"-Sud", varOut[i], youtS_, yscalS_,
                       flag_rms=flag_rms, flag_scat=flag_scat, flag_histo=flag_histo,
                       nbins=nbins, calX=calX);
            
        min_= np.min(y_scale[i]);  max_= np.max(y_scale[i])
        moy_= np.mean(y_scale[i]); std_=np.std(y_scale[i])
        print("%s scalled : %s_decodee(out9), min=%f, max=%f, mean=%f, std=%f"
                    %(strset, varOut[i], min_, max_, moy_, std_));                   

        if not FIGBYPASS and len(im2show) > 0 :
            #print("taille Y_scale:",np(y_scale))
            suptitre="some %s %s_Yscalled decodee %s"%(strset, varOut[i],im2show);

            showsome(y_scale[i][im2show,0,:,:], wmin=wbmin[wk_], wmax=wbmax[wk_],
                     wnorm=wnorm[wk_], cmap=wcmap[wk_], calX0=calX[0][im2show], Xref=VXout_brute[i][im2show,0,:,:],
                     wdif="", wdifr=None, varlib=varOut[i], suptitre=suptitre,savefig=savefig); 
                                      

      
    # R�sultat sur vecteur uv ; hors de la boucle parce qu'on en a besoin que
    # d'une fois, et qu'il faut attendre la fin de la boucle pour que
    # y_scale soit d�cod�; ATTENTION, cela n'est fait que pour les targetout !!!
    resultuv(VXout_brute,y_scale,strset, flag_cosim=flag_cosim, 
             flag_nrj=flag_nrj, flag_ens=flag_ens, flag_corrplex=flag_corrplex,
             flag_div=flag_div, dx=dxRout, dy=dyRout, flag_scat=flag_scat,
             flag_histo=flag_histo, calX=calX, wdif='log');
    #                                       
    if not FIGBYPASS :
        if IS_HUVout and 1 : # Vecteurs UV on H    
            suptitre = "some %s(out) H brute scalled(+uv) %s"%(strset,im2show); #! (H en dure ?)
            ih_ = showhquiv(y_scale, im2show, calX0=calX[0], Yref=VXout_brute, wdif="", 
                            wdifr=None, suptitre=suptitre,savefig=savefig);
#%%
#------------
print("Taille de x_train: ",len(x_train))
if FLGA_RES : # sur APP
    setresult(Mdl, x_train, "APP", VAout_brute, calA, flag_rms=FLGA_RMS,
              flag_scat=FLGA_SCAT, flag_histo=FLGA_HISTO, nbins=NBINS,
              flag_cosim=FLGA_COSIM, flag_nrj=FLGA_NRJ, flag_ens=FLGA_ENS,
              flag_corrplex=FLGA_CORRPLEX, flag_div=FLGA_DIV,
              dx=dxRout, dy=dyRout, wdif='log', im2show=im2showA);
#------------


#------------
print("Taille de x_test: ",len(x_test))
#------------
        
#for i in np.arange(NvarIn):
#    x_test[i] = x_test[i].transpose(0,1,3,2)
#for i in np.arange(NvarOut):
#    y_test[i] = y_test[i].transpose(0,1,3,2)

#%%
if TEST_ON and FLGT_RES : # sur TEST
    setresult(Mdl, x_test, "TEST", VTout_brute, calT, flag_rms=FLGT_RMS,
              flag_scat=FLGT_SCAT, flag_histo=FLGT_HISTO, nbins=NBINS,
              flag_cosim=FLGT_COSIM, flag_nrj=FLGT_NRJ, flag_ens=FLGT_ENS,
              flag_corrplex=FLGT_CORRPLEX, flag_div=FLGT_DIV,
              dx=dxRout, dy=dyRout, wdif='log',  im2show=im2showT,savefig=SAVEFIG);
    #setresult(Mdl, x_test, "TEST", y_test, calT, flag_rms=FLGT_RMS,
    #          flag_scat=FLGT_SCAT, flag_histo=FLGT_HISTO, nbins=NBINS,
    #          flag_cosim=FLGT_COSIM, flag_nrj=FLGT_NRJ, flag_ens=FLGT_ENS,
    #          flag_corrplex=FLGT_CORRPLEX, flag_div=FLGT_DIV,
    #          dx=dxRout, dy=dyRout, wdif='log',  im2show=im2showT);
#------------
plt.show();
#======================================================================
#%%
