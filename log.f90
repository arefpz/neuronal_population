MODULE log
    USE VARS
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE write_log()
    OPEN(UNIT = 10, FILE = "log.txt", ACTION = "WRITE", STATUS = "UNKNOWN")
    WRITE(10, *) "Variable=", "Value"
    WRITE(10, *) "STDP =", STDP
    WRITE(10, *) "STDPEE =", STDPEE
    WRITE(10, *) "STDPEI =", STDPEI
    WRITE(10, *) "STDPIE =", STDPIE
    WRITE(10, *) "STDPII =", STDPII   ! STDP ON/OFF
    WRITE(10, *) "STDPEEH =", STDPEEH 
    WRITE(10, *) "STDPEIH =", STDPEIH
    WRITE(10, *) "STPIEH =", STDPIEH
    WRITE(10, *) "STDPIIH =", STDPIIH ! HEBBIAN STDP
    
    WRITE(10, *) "ompnum =", ompnum
    WRITE(10, *) "refractory =",  refractory
    WRITE(10, *) "Neurons =", Neurons
    WRITE(10, *) "Eneuron=", Eneuron
    WRITE(10, *) "Ineuron=", Ineuron
    WRITE(10, *) "ensemble =", ensemble
    WRITE(10, *) "ptime =", ptime
    WRITE(10, *) "omega = ", omega
    WRITE(10, *) "ampilude=", amp
    WRITE(10, *) "delat_ampilude=", damp
    WRITE(10, *) "IS_Signal=", IS_Signal
    WRITE(10, *) "start_signal=", start_signal
    WRITE(10, *) "end_signal=", end_signal
        
    WRITE(10, *) "sec =", sec
    WRITE(10, *) "InputE =", InputE
    WRITE(10, *) "InputI =", InputI
    WRITE(10, *) "Ein =", Ein
    WRITE(10, *) "Iin =", Iin
    WRITE(10, *) "Esynout =", Esynout
    WRITE(10, *) "EsynE =", EsynE
    WRITE(10, *) "EsynI =", EsynI
    WRITE(10, *) "gfactor =", gfactor
    WRITE(10, *) "gee =", gee
    WRITE(10, *) "gei =", gei
    WRITE(10, *) "gie =", gie
    WRITE(10, *) "gii =", gii
    WRITE(10, *) "dgee =", dgee
    WRITE(10, *) "dgei =", dgei
    WRITE(10, *) "dgie =", dgie
    WRITE(10, *) "dgii =", dgii

    WRITE(10, *) "geemin =", geemin
    WRITE(10, *) "geemax =", geemax
    WRITE(10, *) "geimin =", geimin
    WRITE(10, *) "geimax =", geimax
    WRITE(10, *) "giemin =", giemin
    WRITE(10, *) "giemax =", giemax
    WRITE(10, *) "giimin =", giimin
    WRITE(10, *) "giimax =", giimax
    
    WRITE(10, *) "tau_gaba_r =", tau_gaba_r
    WRITE(10, *) "tau_ampa_r =", tau_ampa_r
    WRITE(10, *) "tau_gaba_d =", tau_gaba_d
    WRITE(10, *) "tau_ampa_d =", tau_ampa_d
    
    WRITE(10, *) "dtau_gaba_r =", dtau_gaba_r
    WRITE(10, *) "dtau_ampa_r =", dtau_ampa_r
    WRITE(10, *) "dtau_gaba_d =", dtau_gaba_d
    WRITE(10, *) "dtau_ampa_d =", dtau_ampa_d

    WRITE(10, *) "Aplus =", Aplus
    WRITE(10, *) "Aminus =", Aminus
    WRITE(10, *) "tauPlus =", tauPlus
    WRITE(10, *) "tauMinus =", tauMinus
    WRITE(10, *) "tauPlusIE =", tauPlusIE
    
    WRITE(10, *) "VrestE =",VrestE
    WRITE(10, *) "dvrestE =", dvrestE
    WRITE(10, *) "TAUME =", TAUME
    WRITE(10, *) "dtaumE =", dtaumE
    WRITE(10, *) "VresetE =", VresetE
    WRITE(10, *) "dvresetE =", dvresetE
    WRITE(10, *) "VTHEE =", VTHEE 
    WRITE(10, *) "dtheE =", dvtheE
    
    WRITE(10, *) "VrestI =", VrestI
    WRITE(10, *) "dvrestI =", dvrestI
    WRITE(10, *) "TAUMI =", TAUMI
    WRITE(10, *) "dtaumI =", dtaumI
    WRITE(10, *) "VresetI =", VresetI
    WRITE(10, *) "dvresetI =", dvresetI
    WRITE(10, *) "VTHEI =", VTHEI
    WRITE(10, *) "dvtheI =", dvtheI
    
    CLOSE(UNIT = 10)
    END SUBROUTINE write_log
END MODULE log
