MODULE vars
    IMPLICIT NONE
    INTEGER, PARAMETER:: DP = selected_real_kind(6, 37)
    ! OPEN MP
    INTEGER:: ompnum = 4
    
    ! TIME
    INTEGER :: cr, c1, c2
    REAL :: systime, time_rate
    
    ! Simulation
    REAL(KIND = DP), PARAMETER  :: sec = 30 ! Simulation time in seconds
    INTEGER :: time
    INTEGER, PARAMETER :: ms = 1000 ! Each sec is ms mili-seconds
    INTEGER :: ensemble= 5 ! Number of trials
    INTEGER :: part = 0, ptime = 10000  ! to write data in several part each ptime steps
    LOGICAL :: PRINT_EACH_PART = .TRUE.
    REAL(KIND = DP), PARAMETER  :: dt = 0.1_DP ! Time step used to integerate the differential equation
    INTEGER :: nsave = 50 ! Number of saved neurons
    
    ! Signal: sinosuidal signal with "omega" Frequency, and "amp * damp" Amplitude, as tACS
    REAL(KIND = DP):: omega = 25.0_DP, amp = 2.0_DP, damp = 0.1_DP
    INTEGER(KIND = 1), PARAMETER :: IS_Signal = 1 ! Is the stimulation presented 1, or not 0
    INTEGER(KIND = DP), PARAMETER :: start_signal=5, end_signal=40 ! Starting and ending time point for stimulation
    
    ! Neuron Parameters
    INTEGER :: Neurons = 10000 ! Number of neurons
    INTEGER :: Eneuron, Ineuron
    INTEGER, PARAMETER :: refractory = INT(2/dt) ! Refractory period in time step based on dt
    REAL, PARAMETER  :: NeRate = 0.8_DP ! Excitatory to Inhibitory neurons Ratio (E out of N). 0.8 is 4:1 E:I
    REAL(KIND = DP), PARAMETER  :: Esynout = 0.0_DP ! Reversal potential of E neurons recieving Poisonnian spike train.
    REAL(KIND = DP), PARAMETER  :: EsynE = 0.0_DP, EsynI = -85.0_DP ! Reversal Potential for E and I neurons
    
    REAL(KIND = DP), PARAMETER  :: VrestE  = -60.0_DP,  dvrestE = 0.0_DP ! Resting potential for LIF neurons and its sigma
    REAL(KIND = DP), PARAMETER  :: TAUME   = 10.0_DP,    dtaumE = 3_DP
    REAL(KIND = DP), PARAMETER  :: VresetE = -60.0_DP, dvresetE = 0.0_DP
    REAL(KIND = DP), PARAMETER  :: VTHEE   = -54.0_DP,   dvtheE  = 0.0_DP
    
    REAL(KIND = DP), PARAMETER  :: VrestI  = -60.0_DP,  dvrestI = 0.0_DP
    REAL(KIND = DP), PARAMETER  :: TAUMI   = 10.0_DP,    dtaumI = 3.0_DP
    REAL(KIND = DP), PARAMETER  :: VresetI = -60.0_DP, dvresetI = 0.0_DP
    REAL(KIND = DP), PARAMETER  :: VTHEI   = -54.0_DP,   dvtheI = 0.0_DP

    REAL(KIND = DP):: InputE = 5.5_DP, InputI = 5.5_DP, ISD = 1_DP
    
    !! NETWORK PARAMETERS
    ! Probability of Random Connectivity 
    REAL(KIND = DP), PARAMETER  :: PEE = 0.1_DP, PEI = 0.1_DP, PIE = 0.1_DP, PII = 0.1_DP
    ! Axonal Delay in ms
    REAL(KIND = DP), PARAMETER  :: ax_dee = 0.0_DP, ax_dei = 0.0_DP, ax_die = 0.0_DP, ax_dii = 0.0_DP
    REAL(KIND = DP), PARAMETER  :: ax_ddee = 1.0_DP, ax_ddei = 1.0_DP, ax_ddie = 0.5_DP, ax_ddii = 0.5_DP
    
    !! STNAPSES PARAMETERS
    ! Synaptic Weights
    REAL(KIND = DP), PARAMETER  :: gfactor = 0.0005 ! for N = 10000, gfactor is 0.2
    REAL(KIND = DP), PARAMETER  :: gee = gfactor * 0.1_DP
    REAL(KIND = DP), PARAMETER  :: gei = gfactor * 0.1_DP
    REAL(KIND = DP), PARAMETER  :: gie = gfactor * 0.5_DP
    REAL(KIND = DP), PARAMETER  :: gii = gfactor * 0.5_DP
    REAL(KIND = DP), PARAMETER  :: dgee = 0.1*gee, dgei = 0.1*gei, dgie = 0.1*gie, dgii = 0.1*gii
    ! Synaptic Efficacy Parameters (Double exponential function)
    REAL(KIND = DP), PARAMETER  :: tau_gaba_r = 0.5_DP, dtau_gaba_r = 0.1_DP
    REAL(KIND = DP), PARAMETER  :: tau_ampa_r = 0.5_DP, dtau_ampa_r = 0.1_DP
    REAL(KIND = DP), PARAMETER  :: tau_gaba_d = 5.0_DP, dtau_gaba_d = 0.5_DP
    REAL(KIND = DP), PARAMETER  :: tau_ampa_d = 3.0_DP, dtau_ampa_d = 0.5_DP
    
    ! STDP PARAMETERS
    LOGICAL :: STDP = .TRUE.
    LOGICAL :: STDPEE = .TRUE., STDPEI = .TRUE., STDPIE = .TRUE., STDPII = .FALSE.   ! STDP ON/OFF
    LOGICAL :: STDPEEH = .TRUE., STDPEIH = .TRUE., STDPIEH = .TRUE., STDPIIH = .FALSE. ! HEBBIAN STDP
    
    REAL(KIND = DP), PARAMETER  :: geemin = 0.01_DP * gee, geemax = 2.0_DP * gee
    REAL(KIND = DP), PARAMETER  :: geimin = 0.01_DP * gei, geimax = 2.0_DP * gei
    REAL(KIND = DP), PARAMETER  :: giemin = 0.01_DP * gie, giemax = 2.0_DP * gie
    REAL(KIND = DP), PARAMETER  :: giimin = 0.01_DP * gii, giimax = 2.0_DP * gii
    
    REAL(KIND = DP), PARAMETER :: Aplus = 0.4_DP * gee, Aminus = 0.2_DP * gee
    REAL(KIND = DP), PARAMETER :: tauPlus = 10.0_DP, tauMinus = 10.0_DP, tauPlusIE = 30
    
    ! VARIABLES
    integer, allocatable :: seed(:)
    INTEGER :: nseed
    REAL(KIND = DP):: PI = 4*ATAN(1.0_DP)
    
    INTEGER :: ii, jj, ll, kk, mm, ierr, max_delay_condition, pti, ngei, ngie
    INTEGER :: ens
    INTEGER :: dI0
    INTEGER :: maxCon
    INTEGER, DIMENSION(:), ALLOCATABLE :: st1, st2
    INTEGER, DIMENSION(:), SAVE, ALLOCATABLE :: ABnum, EABnum, IABnum  ! sigma  ! XA, XD, 
    INTEGER(KIND = 1), DIMENSION(:), ALLOCATABLE :: fired
    INTEGER, DIMENSION(:), ALLOCATABLE :: tref
    INTEGER(KIND = 1), DIMENSION(:,:), ALLOCATABLE :: A
    INTEGER, DIMENSION(:,:), SAVE, ALLOCATABLE :: AB, ABdelay_ax, ABdelay_den
    INTEGER(KIND = 1), DIMENSION(:,:), ALLOCATABLE :: rho
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: delay_ax
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: delay_den
    INTEGER(KIND = 4), DIMENSION(:,:), ALLOCATABLE :: selected_synapses
    INTEGER(KIND = 4), DIMENSION(:), ALLOCATABLE :: selected_neurons
    INTEGER :: delay_max
    INTEGER :: SPt = 0  ! If the first part has pased 1 or not 0. Initially it should be 0
    INTEGER :: sec1  ! 1 second
    
    REAL :: Ein = 0, Iin = 0
    REAL :: start1, start2, finish1, finish2
    REAL(KIND = DP) :: AuxSum
    REAL(KIND = DP) :: deltaT, dweight
    REAL(KIND = DP) :: tempJSumE, tempJSumI
    REAL(KIND = DP) :: xv, x, xt
    REAL(KIND = DP) :: xsyn, xnoise, S1, S2!, S3
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: vth, vreset, vrest
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: S, Isyn
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: v
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: taum
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: taurise, taudecay 
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: Esyn, EsynSign
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: JE0
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: JI0
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: synapse
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: SS, SSx
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Vsave
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: ABsynapse
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: vrnd
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: xIsyn
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: Iext
    
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: VEneu, VIneu  
    REAL(KIND = DP), ALLOCATABLE, DIMENSION(:,:) :: Isyn_save, stdp_save, IsynEI 
    REAL(KIND = DP), ALLOCATABLE, DIMENSION(:) :: Isignal
    REAL(KIND = DP), ALLOCATABLE, DIMENSION(:) :: I0
    REAL(KIND = 4), ALLOCATABLE, DIMENSION(:) :: gsumee, gsumei, gsumie
    

    REAL(KIND = DP), ALLOCATABLE, DIMENSION(:):: Inoise
    ! REAL(KIND = DP):: xsum1, xsum2
    CHARACTER(len = 1024) :: filename1
    INTEGER(KIND = DP), ALLOCATABLE, DIMENSION(:) :: Xpoissonsp
    
END MODULE vars
