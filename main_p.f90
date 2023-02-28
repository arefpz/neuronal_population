!! Compiled on 21 Feb. 2021 by gfortran 10.2.0 on Ubuntu 20.10 as following lines
! Non-Parallel
!! gfortran -m64 -O3 vars.f90 log.f90 writer.f90 random.f90 modu.f90 main_p.f90 -o a.out
! Parallel with openmp
!! gfortran -m64 -O3 -fopenmp vars.f90 log.f90 writer.f90 random.f90 modu.f90 main_p.f90 -o a.out

!! echo <$var1> <$var2> | ./a.out ! No need to <>. 

!! Written by Aref Pariz & Farhad Daei, Email: pariz.aref@gmail.com, farhad.daei@gmail.com
!! The full version will be available on developers' github page ASAP
PROGRAM LIFpop

    USE iso_fortran_env
    USE VARS
    USE MODU
    USE omp_lib
    USE WRITER    
    USE LOG
    USE random
    IMPLICIT NONE
    integer:: endIndex, indexTemp
    real(kind = DP):: sumReduction
    real(kind = DP), save, dimension(:), allocatable:: temp1, temp2
    !$omp threadprivate (temp1, temp2)
    time = INT(1000*sec/dt)
    
    Eneuron = INT(NeRate*neurons) 
    Ineuron = INT(neurons-Eneuron)
    maxCon = INT(MAXVAL([Eneuron*PEE+Ineuron*PIE, Eneuron*PEI+Ineuron*PII]))
    time = time+mod(8-mod(time, 8), 8)
    print*, "Max Con :", maxCon
    print*, "Time :", time
    print*, "E, I :", Eneuron, Ineuron
    CALL system_clock(count_rate=cr)
    time_rate = REAL(cr)
    CALL SYSTEM_CLOCK(c1)

    IF (mod(neurons, 8) /= 0) THEN
        write(0, *) "Error:: Number of Neurons must be divisiable by 8 due to performance "
        STOP
    END IF
    !dir$ assume (mod(time, 8).eq.0)
    !dir$ assume (mod(neurons, 8).eq.0)
    !$OMP PARALLEL &
    !$OMP NUM_THREADS(ompnum)
    ALLOCATE(temp1(neurons), temp2(neurons))
    !$OMP END PARALLEL
    ALLOCATE(synapse(neurons, neurons), & ! Synaptic weight matrix
            fired(neurons), & ! Fired neruons =1, else=0 
            tref(neurons), & ! Refractory time steps
            taudecay(neurons), & ! synaptic current efficacy decay time
            taurise(neurons), & ! synaptic current efficacy rise time
            !C(neurons), & ! Capacitance. Not used Here
            I0(neurons), & ! Constant input current
            vth(neurons), & ! Threshold for spike generation
            A(neurons, neurons), & ! Connectivity matrix
            delay_ax(neurons, neurons), & ! Axonal delay matrix
            Isyn(neurons ), & ! Synaptic current 
            Isignal(time), & ! Signal, tACS
            Inoise(neurons), & ! Noise input
            xIsyn(maxCon), & ! Auxiliary synaptic input for each neuron
            v(neurons ), & ! Neurons' membrane potential
            taum(neurons ), & ! membrane time Constant
            vreset(neurons ), & ! Reseting potential
            vrest(neurons ), & ! Resting potential
            Esyn(neurons), & ! Synapse, Reveral potential
            EsynSign(neurons), & ! Sign of synaptic reversal potential, positive -> excitaotry, negative-> inhibitory
            EABnum(neurons), & ! Excitatory pre-synaptic neurons index
            IABnum(neurons), & ! inhibitory pre-synaptic neurons index
            S(neurons ), & ! Synaptic efficacy
            gsumee(ptime), & ! Sum of E->E synaptic weights
            gsumei(ptime), & ! Sum of E->I synaptic weights
            gsumie(ptime), & ! Sum of I->E synaptic weights
            st1(neurons ), & ! last spike effect
            st2(neurons ), & ! one before last spike effect
            ABnum(neurons ), & ! 
            JE0(neurons ), &
            JI0(neurons ), &
            Iext(neurons), &
            selected_neurons(nsave), &
            selected_synapses(nsave,2), &
            Isyn_Save(nsave, ptime), &
            stdp_Save(nsave, ptime), &
            VEneu(ptime), VIneu(ptime),&
            Vsave(nsave, ptime), &
            rho(neurons, ptime), &
            AB(maxCon, neurons), &
            ABsynapse(maxCon, neurons), &
            ABdelay_ax(maxCon, neurons), &
            ABdelay_den(maxCon, neurons), &
            vrnd(neurons), &
            STAT = ierr)
    IF (ierr /= 0) THEN
            PRINT*, "Error:: in ALLOCATION, check the DIMENSION of matrices"
            STOP
    END IF
    
    PRINT*, "Allocation succeed"
    CALL cpu_time(start1)
    READ(*,*) omega, amp
    OPEN(UNIT = 10, FILE = "variable.txt", ACTION = "WRITE", STATUS = "UNKNOWN")
    WRITE(10, *) omega, amp, STDPEE, STDPEEH, STDPEI, STDPEIH, STDPIE, STDPIEH
    CLOSE(UNIT = 10)
    IF (STDP .EQV. .TRUE.) THEN
        PRINT*, "STDP= ", STDP
        IF (STDPEE .EQV. .TRUE.) THEN
            PRINT*, "STDP E2E is Hebbian= ", STDPEEH
        END IF
        IF (STDPEI .EQV. .TRUE.) THEN
            PRINT*, "STDP E2I is Hebbian= ", STDPEIH
        END IF
        IF (STDPIE .EQV. .TRUE.) THEN
            PRINT*, "STDP I2E is Hebbian= ", STDPIEH
        END IF
        IF (STDPII .EQV. .TRUE.) THEN 
            PRINT*, "STDP I2I is Hebbian= ", STDPII
        END IF

    END IF
    IF (IS_Signal == 1) THEN
        PRINT*, "Omega= ", omega
        PRINT*, "Amplitude= ", amp * damp, "mV"
    ELSE
        PRINT*, "No Signal"
    END IF
    CALL write_log()

!     CALL RANDOM_SEED(SIZE = nseed)
!     ALLOCATE(seed(nseed))
!     seed = [5465464,43489,98786,908324,142876,21095,232098,3010029]
!     CALL RANDOM_SEED(put=seed)

    !!-----------------Periodic Signal-------------START----------------
    Isignal(:) = 0
    IF (IS_Signal == 1) THEN
        DO ii = 0, time-1
            Isignal(ii+1) = damp * amp * SIN( omega*2 * PI*ii*dt/ms )
        END DO
        sec1 = INT(1*ms/dt)
        Isignal(0*sec1+1:start_signal*sec1) = 0
        Isignal(end_signal*sec1+1:time) = 0
        PRINT*, "Isignal got Created"
        CALL write_1DR(Isignal,"Isignal")

    END IF
    !!-----------------Periodic Signal-------------END------------------
    
    !!-%%%%%%%%%%%%%%%%%%-MAIN PART OF PROGRAM STARS-%%%%%%%%%%%%%%%%%%%%%%%-
    DO ens =1, ensemble
        pti = 0
        part = 0
        
        selected_neurons = intrnd(1, Eneuron, nsave)
        selected_synapses(:,1) = intrnd(1, maxCon, nsave)  ! PRE
        selected_synapses(:,2) = intrnd(1, Eneuron, nsave) ! POST
        CALL write_1DI(selected_neurons,"selcN")
        CALL write_2DI(selected_synapses,"selcS")
        
    
        CALL cpu_time(start2)
!!------------------INITIALIZATION of MATRICES---------------------
        A = 0
        delay_ax = 0
        synapse = 0.0_DP
        max_delay_condition = 1
        VIneu = 0.0_DP
        VEneu = 0.0_DP
        Vsave = 0.0_DP
        stdp_Save = 0.0_DP
        Isyn_Save = 0.0_DP
        gsumee = 0.0_DP
        gsumei = 0.0_DP
        gsumie = 0.0_DP

!!-----------------GENERATING CONCCETION MATRIX START---------------

        A = creatA_fixed_number(neurons, Eneuron, Ineuron, Pee, Pei, Pie, Pii)
        DO ii = 1, neurons
            A(ii, ii) = 0
        END DO
        PRINT*, "Connectivity Matrix got created"
!!-----------------GENERATING CONNECTION MATRIX END-----------------

!!------------------GENERATING DELAY SECTION----START---------------

        delay_ax = create_delay(ax_dee, ax_dei, ax_die, ax_dii, &
                                ax_ddee, ax_ddei, ax_ddie, ax_ddii,&
                                neurons, Eneuron, Ineuron, dt, A)
        DO ii = 1, neurons
            DO jj = 1, neurons
                IF (A(ii, jj) .EQ. 0) THEN
                    delay_ax(ii, jj) = 0
                    ! delay_den(ii, jj) = 0
                END IF
            END DO
        END DO
        delay_max = maxval(delay_ax)
        ALLOCATE(SS(neurons, ptime+delay_max), SSx(neurons, delay_max))
        PRINT*, size(SS,1), size(SS,2)
        print*, "Axonal Delay Matrix got created"
!!------------------GENERATING DELAY SECTION----END-----------------

!!------------------SYNAPTIC WEIGHT SECTION----START----------------
        synapse = synapses(gee, gei, gie, gii, dgee, dgei, dgie, dgii, neurons, Eneuron, Ineuron, A)
        DO ii = 1, neurons
            DO jj = 1, neurons
                IF (A(ii, jj) .EQ. 0) THEN
                    synapse(ii, jj)=0
                END IF
            END DO
        END DO
        PRINT*, "Synaptic Weight Matrix got created"
!         CALL write_2DR(synapse, 'synapse')

!!------------------SYNAPTIC WEIGHT SECTION---------END------------------

!!------------------GENERATION INOISE--------------START----------------
!         DO ii = 1, time
!            DO ll = 1, neurons
!                 xnoise = random_normal()
!                 Inoise(ll, ii) = sqrt(ISD) * xnoise / sqrt(dt)
!            END DO
!         END DO
!!------------------GENERATION INOISE-------------END------------------        
!!-----------------INITIALIZATION------------------START----------     
        rho = 0
        tref = 0
        SS = 0.0_DP
        SSx = 0.0_DP
        S1 = 0.0_DP
        S2 = 0.0_DP
        fired = 0
        
        !S3 = 0
        
        st1 = -100000
        st2 = -100000
        ! st3 = -100000
        dI0 = 0
        I0(1:Eneuron) = InputE+Ein*0.025_DP
        I0(Eneuron+1:Neurons) = InputI+Iin*0.025_DP
    
        Esyn(1:Eneuron) = EsynE
        Esyn(Eneuron+1:Neurons) = EsynI
        
        EsynSign(1:Eneuron) = 1.0_DP
        EsynSign(Eneuron+1:Neurons) = -1.0_DP
        
        DO ii = 1, Eneuron
            xnoise = random_normal()
            taudecay(ii) = tau_ampa_d + SQRT(dtau_ampa_d) * xnoise;  ! AMPA  DECAY TIME
            call random_number(xnoise)
            taurise(ii) = tau_ampa_r;                   ! AMPA  RISE  TIME
            
        END DO
        DO ii = Eneuron + 1, Neurons
            xnoise = random_normal()
            taudecay(ii) = tau_gaba_d + SQRT(dtau_gaba_d) * xnoise;   ! GABAa DECAY TIME
            call random_number(xnoise)
            taurise(ii) = tau_gaba_r;    ! GABAa RISE  TIME
        END DO

        call random_number(vrnd)
        v = -60.0_DP + vrnd
        Isyn = 0.0_DP
        Iext = 0.0_DP
        DO ii = 1, Eneuron
            xnoise = random_normal()
            vth(ii) = VTHEE + SQRT(dvtheE) * xnoise
            xnoise = random_normal()
            vreset(ii) = VresetE + SQRT(dvresetE) * xnoise
            xnoise = random_normal()
            taum(ii) = TAUME + SQRT(dtaumE) * xnoise
            xnoise = random_normal()
            vrest(ii) = VrestE + SQRT(dvrestE) * xnoise
        END DO
        DO ii = Eneuron+1, Neurons
            xnoise = random_normal()
            vth(ii) = VTHEI + SQRT(dvtheI) * xnoise
            xnoise = random_normal()
            vreset(ii) = VresetI + SQRT(dvresetI) * xnoise
            xnoise = random_normal()
            taum(ii) = TAUMI + SQRT(dtaumI) * xnoise
            xnoise = random_normal()
            vrest(ii) = VrestI + SQRT(dvrestI) * xnoise
        END DO
        AB = 0
        ABnum = 0
        EABnum = 0
        IABnum = 0
        ABdelay_ax = 0
        ABdelay_den = 0
        ABsynapse = 0.0_DP
        stdp_Save = 0.0_DP
        JE0 = 0.0_DP
        JI0 = 0.0_DP
        ! Each cloumn is the post synaptic neruon, and each row is pre synaptic neuron.
        ! i.e. ABsynapse(10, 5) means the neuron 5 recives synaptic input from neuron AB(10, 5)
        DO jj = 1, Neurons
            kk = 0
            DO ii = 1, Neurons
                IF (A(ii, jj) .EQ. 1) THEN
                    kk = kk+1
                    AB(kk, jj) = ii !! "kk'th" neighbor of neuron number "jj", is the neuron number "ii"
                    ABdelay_ax(kk, jj) = delay_ax(ii, jj) !! the axonal time delay between neuron number "jj",
                                                            ! and it's kk'th" neighbor (which is the neuron number "ii")
                    ! ABdelay_den(kk, jj) = delay_den(ii, jj) !! the dendritic time delay between neuron number "jj"
                                                            !, and it's kk'th" neighbor (which is the neuron number "ii")
                    ABsynapse(kk, jj) = synapse(ii, jj)  !! same shit as the line above, but for synapse; -)
                    IF (ii .LE. Eneuron) THEN
                        EABnum(jj) = kk
                    ELSE
                        IABnum(jj) = kk - EABnum(jj)
                    END IF
                END IF
            END DO
            ABnum(jj) = kk
        END DO

        DO jj=1, Neurons
            JE0(jj) = SUM(ABsynapse(1:EABnum(jj),jj))
            JI0(jj) = SUM(ABsynapse(EABnum(jj)+1:EABnum(jj)+IABnum(jj),jj))
        END DO
        CALL write_2DI(AB,'ABval')
        CALL write_2DR(ABsynapse,'ABsynapse0')
        CALL write_1DR(taum,'taum')
        CALL write_1DR(vth,'vthreshold')
        CALL write_1DR(vrest,'vrest')
        CALL write_1DR(vreset,'vreset')


        !! ----------LOOP OVER TIME AND NEURONS--------START----------------
        PRINT*, "Main part is started"
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP SHARED(S, Isyn, time, maxCon, max_delay_condition, VEneu, VIneu) &
        !$OMP SHARED(selected_neurons, selected_synapses, Isyn_Save, Vsave, nsave, stdp_Save) &
        !$OMP SHARED(rho, tref, st1, st2, ptime) &
        !$OMP SHARED(Inoise, vth, taurise, taudecay, vreset, vrest, taum, I0) &
        !$OMP SHARED(ABnum, ABdelay_ax, ABdelay_den, ABsynapse, EABnum, IABnum, SS, SSx, Esyn) &
        !$OMP SHARED(v, AB, neurons, fired, Eneuron, delay_max, spt, part, EsynSign, STDP) &
        !$OMP SHARED(STDPEE, STDPEI, STDPIE, STDPEEH, STDPEIH, STDPIEH) &
        !$OMP SHARED(PRINT_EACH_PART), &
        !$OMP SHARED(Isignal, ISD), &
        !$OMP SHARED(gsumee, gsumei, gsumie, AuxSum), &
        !$OMP PRIVATE(xIsyn, deltaT, dweight, endIndex, indexTemp, sumReduction)&
        !$OMP PRIVATE(S1, S2, kk, jj, ll, xt, pti, ii, mm)&
        !$OMP PRIVATE(tempJSumE, tempJSumI)&
        !$OMP NUM_THREADS(ompnum)
        pti = 0
        DO ii = 1, time
            !$OMP SINGLE
            DO ll = 1, neurons
                Inoise(ll) = sqrt(ISD) * random_normal() / sqrt(dt)
            END DO
            !$OMP END SINGLE
            !$OMP DO
            DO ll = 1, Neurons
                IF (ii .GT. tref(ll)) THEN
                    ! Euler-Maruyama Method
                    v(ll) = v(ll) + dt*funv(v(ll), vrest(ll), Isyn(ll), Inoise(ll) + &
                    I0(ll) + Isignal(ii), taum(ll))

                    IF (v(ll) .GT. vth(ll)) THEN
                        st1(ll) = st2(ll)  ! Saveing the previous spike time
                        st2(ll) = ii   ! saveing the current spike time
                        fired(ll) = 1  ! the neuron spiked at ii
                        tref(ll) = ii + refractory  ! the membrane should stay fixed for refractory period
                        v(ll) = vreset(ll)  ! Keeping the v at vreset for refractory period
                    ELSE
                        fired(ll) = 0
                        
                    END IF
                ELSE
                    v(ll) = vreset(ll)  ! The ii should be more that tref to be able to evolve again
                    fired(ll) = 0 ! 
                END IF
                ! Below is a normalizing factor for synapse efficacy
                xt = (taurise(ll)/taudecay(ll))**(taurise(ll)/(taudecay(ll)-taurise(ll))) &
                -(taurise(ll)/taudecay(ll))**(taudecay(ll)/(taudecay(ll)-taurise(ll)))
                ! spike before last spike
                S1 = (EXP(-dt*(ii-st1(ll))/taurise(ll))-EXP(-dt*(ii-st1(ll))/taudecay(ll))) / xt 
                ! The last spike
                S2 = (EXP(-dt*(ii-st2(ll))/taurise(ll))-EXP(-dt*(ii-st2(ll))/taudecay(ll))) / xt
                
                S(ll) =  S1 + S2
                
            END DO 
            !$OMP END DO
            pti = pti+1  ! counter that resets every ptime
            !$OMP SINGLE
            rho(:,pti) = fired  ! Spike train
            SS(:, pti + (SPt * delay_max)) = S  ! synaptic efficacy of outgoing synapses from pre neruons
            VEneu(pti) = SUM(v(1:Eneuron))
            VIneu(pti) = SUM(v(Eneuron+1:neurons))
            !$OMP END SINGLE
            
            SELECT CASE (max_delay_condition)  ! if the delay less than max_delay
                CASE(1)
                    DO jj = 1, neurons  ! Post
                        xIsyn = 0
                        DO kk = 1, ABnum(jj)  ! Pre
                            IF ((ii .GT. ABdelay_ax(kk, jj))) THEN
                                xIsyn(kk) = SS(AB(kk, jj), ii-ABdelay_ax(kk, jj) ) &
                                 * ABsynapse(kk, jj) * (v(jj) - Esyn(AB(kk, jj)))
                            ELSE
                                xIsyn(kk) = 0
                            END IF
                        END DO
                        Isyn(jj) = SUM(xIsyn)
                    END DO
                CASE(2)
                    IF (STDP .EQV. .TRUE.) THEN
                    ! ======================== STDP =========START==========
                        ! STDP among E - > E neurons
                        IF (STDPEE .EQV. .TRUE.) THEN
                            !$OMP DO
                            DO jj = 1, Eneuron  ! POST
                                DO kk = 1, EABnum(jj)  ! PRE
                                    IF (ii .EQ. (st2(AB(kk, jj)) + ABdelay_ax(kk, jj)))  THEN  ! Pre neuron spikes arrives on Synapse
                                        deltaT = dt * (st2(jj) - ii)
                                    ELSEIF (fired(jj) .EQ. 1) THEN  ! Post neurons spikes
                                        IF (ii .GE. st2(AB(kk,jj)) + ABdelay_ax(kk,jj)) THEN
                                            deltaT = dt * (ii-st2(AB(kk,jj)) - ABdelay_ax(kk, jj))
                                        ELSEIF (ii .LT. st2(AB(kk,jj)) + ABdelay_ax(kk,jj)) THEN
                                            deltaT = dt * (ii-st1(AB(kk,jj)) - ABdelay_ax(kk, jj))
                                        END IF
                                    ELSE
                                        CYCLE
                                    END IF

                                    dweight = 0
                                    IF (STDPEEH .EQV. .TRUE.) THEN
                                        IF (deltaT .GE. 0) THEN
                                            dweight =  Aplus * (1 - ABsynapse(kk,jj)/geemax) * EXP(-deltaT/tauPlus)
                                        ELSE
                                            dweight = - Aminus * (ABsynapse(kk,jj)/gee) * EXP(deltaT/tauMinus)
                                        END IF
                                    ELSE
                                        IF (deltaT .GE. 0) THEN
                                            dweight = - Aminus * (ABsynapse(kk,jj)/gee) * EXP(-deltaT/tauMinus)
                                        ELSE
                                            dweight = + Aplus * (1 - ABsynapse(kk,jj)/geemax) * EXP(deltaT/tauPlus)
                                        END IF
                                    END IF
                                    ABsynapse(kk, jj) = ABsynapse(kk, jj)  + dweight
                                    IF ((ABsynapse(kk, jj)) .GE. geemax) THEN
                                        ABsynapse(kk, jj) = geemax
                                    ELSEIF ((ABsynapse(kk, jj)) .LE. geemin) THEN
                                        ABsynapse(kk, jj) = geemin
                                    END IF  
                                END DO
                            END DO
                            !$OMP END DO 
                        END IF

                        ! STDP among E - > I neurons
                        IF (STDPEI .EQV. .TRUE.) THEN
                            !$OMP DO
                            DO jj = Eneuron + 1, neurons  ! POST
                                DO kk = 1, EABnum(jj)  ! PRE
                                    IF (ii .EQ. (st2(AB(kk, jj)) + ABdelay_ax(kk, jj)))  THEN  ! Pre neuron spikes arrives on Synapse
                                        deltaT = dt * (st2(jj) - ii)
                                    ELSEIF (fired(jj) .EQ. 1) THEN  ! Post neurons spikes
                                        IF (ii .GE. st2(AB(kk,jj)) + ABdelay_ax(kk,jj)) THEN
                                            deltaT = dt * (ii-st2(AB(kk,jj)) - ABdelay_ax(kk, jj))
                                        ELSEIF (ii .LT. st2(AB(kk,jj)) + ABdelay_ax(kk,jj)) THEN
                                            deltaT = dt * (ii-st1(AB(kk,jj)) - ABdelay_ax(kk, jj))
                                        END IF
                                    ELSE
                                        CYCLE
                                    END IF
                                    dweight = 0
                                    IF (STDPEIH .EQV. .TRUE.) THEN
                                        IF (deltaT .GE. 0) THEN
                                            dweight =  Aplus * (1 - ABsynapse(kk,jj)/geimax) * EXP(-deltaT/tauPlus)
                                        ELSE
                                            dweight = - Aminus * (ABsynapse(kk,jj)/gei) * EXP(deltaT/tauMinus)
                                        END IF
                                    ELSE
                                        IF (deltaT .GE. 0) THEN
                                            dweight = - Aminus * (ABsynapse(kk,jj)/gei) * EXP(-deltaT/tauMinus)
                                        ELSE
                                            dweight =   Aplus * (1 - ABsynapse(kk,jj)/geimax) * EXP(deltaT/tauPlus)
                                        END IF
                                    END IF
                                    ABsynapse(kk, jj) = ABsynapse(kk, jj) + dweight
                                    IF ((ABsynapse(kk, jj)) .GE. geimax) THEN
                                        ABsynapse(kk, jj) = geimax
                                    ELSEIF ((ABsynapse(kk, jj)) .LT. geimin) THEN
                                        ABsynapse(kk, jj) = geimin
                                    END IF
                                END DO
                            END DO
                            !$OMP END DO
                        END IF

                        ! STDP among I -> E neurons
                        IF (STDPIE .EQV. .TRUE.) THEN
                            !$OMP DO
                            DO jj = 1, Eneuron  ! POST
                                DO kk = EABnum(jj) + 1, ABnum(jj)  ! PRE
                                    IF (ii .EQ. (st2(AB(kk, jj)) + ABdelay_ax(kk, jj)))  THEN  ! Pre neuron spikes arrives on Synapse
                                        deltaT = dt * (st2(jj) - ii)
                                    ELSEIF (fired(jj) .EQ. 1) THEN  ! Post neurons spikes
                                        IF (ii .GE. st2(AB(kk,jj)) + ABdelay_ax(kk,jj)) THEN
                                            deltaT = dt * (ii-st2(AB(kk,jj)) - ABdelay_ax(kk, jj))
                                        ELSEIF (ii .LT. st2(AB(kk,jj)) + ABdelay_ax(kk,jj)) THEN
                                            deltaT = dt * (ii-st1(AB(kk,jj)) - ABdelay_ax(kk, jj))
                                        END IF
                                    ELSE
                                        CYCLE
                                    END IF
                                    dweight = 0
                                    IF (STDPEIH .EQV. .TRUE.) THEN
                                        IF (deltaT .GE. 0) THEN
                                            dweight = Aplus * (1 - ABsynapse(kk,jj)/giemax) * EXP(-deltaT/tauPlus)
                                        ELSE
                                            dweight = - Aminus * (ABsynapse(kk,jj)/gie) * EXP(deltaT/tauMinus)
                                        END IF
                                    ELSE
                                        IF (deltaT .GE. 0) THEN
                                            dweight = - Aminus * (ABsynapse(kk,jj)/gie) * EXP(-deltaT/tauMinus)
                                        ELSE
                                            dweight = Aplus * (1 - ABsynapse(kk,jj)/giemax) * EXP(deltaT/tauPlus)
                                        END IF
                                    END IF
                                    ABsynapse(kk, jj) = ABsynapse(kk, jj) + dweight
                                    IF (ABsynapse(kk, jj) .GE. giemax) THEN
                                        ABsynapse(kk, jj) = giemax
                                    ELSEIF (ABsynapse(kk, jj) .LE. giemin) THEN
                                        ABsynapse(kk, jj) = giemin
                                    END IF
                                END DO
                            END DO
                            !$OMP END DO
                        END IF
                        
                        ! STDP among I -> I neurons
                        ! Kept fixed
                    END IF
                 
                    !$OMP DO
                    DO jj = 1, neurons  ! POST
                        endIndex = ABnum(jj)
                        temp2(1:endIndex) = Esyn(AB(1:endIndex, jj))
                        DO kk = 1, endIndex  ! PRE
                            indexTemp = pti + (SPt * delay_max) - ABdelay_ax(kk, jj)
                            temp1(kk) = SS(AB(kk, jj), indexTemp)
                        END DO
                        sumReduction = 0.0_DP
                        !$OMP SIMD SIMDLEN(4) ALIGNED(temp1, temp2:64) REDUCTION(+:sumReduction)
                        DO kk = 1, endIndex
                            sumReduction = sumReduction + temp1(kk) * ABsynapse(kk, jj) * (v(jj) - temp2(kk))
                        END DO
                        !$OMP END SIMD
                        Isyn(jj) = sumReduction
                    END DO
                    !$OMP END DO
            END SELECT
            

            !$OMP SINGLE
            AuxSum = 0;
            DO jj=1, Eneuron
                    AuxSum = AuxSum + SUM(ABsynapse(1:EABnum(jj),jj))
            END DO
            gsumee(pti)=AuxSum
            
            AuxSum = 0;
            DO jj=Eneuron + 1, neurons
                    AuxSum = AuxSum + SUM(ABsynapse(1:EABnum(jj),jj))
            END DO
            gsumei(pti)=AuxSum
            
            AuxSum = 0;
            DO jj=1, Eneuron
                    AuxSum = AuxSum + SUM(ABsynapse(EABnum(jj)+1:ABnum(jj),jj))
            END DO
            gsumie(pti)=AuxSum

            DO mm = 1, nsave
                    Isyn_save(mm, pti) = Isyn(selected_neurons(mm))
                    Vsave(mm, pti) = v(selected_neurons(mm))
                    stdp_save(mm, pti) = ABsynapse(selected_synapses(mm,1),selected_synapses(mm,2))
            END DO
            !$OMP END SINGLE
            
            IF (pti .eq. ptime) THEN
                !$OMP SINGLE
                part = part+1
                CALL write_rho_part()
                IF ((part .LE. 10) .OR. (part .GE. 230)) THEN
                    CALL write_2DR(ABsynapse,'ABSynapse')
                END IF
                IF (PRINT_EACH_PART .EQV. .TRUE.) THEN
                    CALL write_2DR(ABsynapse,'ABsynapse')
                   
                    CALL write_2DR(VSave,'Vsave')
                    CALL write_2DR(stdp_save,'stdpc')
                    CALL write_2DR(Isyn_save,'IsynS')
                    
                    CALL write_1DR(VEneu,"VEneu")
                    CALL write_1DR(VIneu,"VIneu")
                    
                    CALL write_1DR(gsumee,"Gsumee")
                    CALL write_1DR(gsumei,"Gsumei")
                    CALL write_1DR(gsumie,"Gsumie")
                    
                    gsumee = 0.0_DP
                    gsumei = 0.0_DP
                    gsumie = 0.0_DP
                    Vsave=0.0_DP
                    stdp_Save=0.0_DP
                    Isyn_Save=0.0_DP
                    VIneu=0.0_DP
                    VEneu=0.0_DP
                END IF
                rho = 0
                SSx = SS(:, ptime + (SPt*delay_max)-delay_max+1: ptime + SPt*delay_max)
                SS = 0
                SS(:, 1: delay_max) = SSx
                !$OMP END SINGLE
                pti = 0
                SPt = 1
            END IF
            
        END DO
        !$OMP END PARALLEL

        DEALLOCATE(SSx, SS)
        CALL cpu_time(finish2)    
        PRINT*, 'Each ensemble time is = ', finish2-start2, ' seconds.'
        !!-----------LOOP OVER TIME AND NEURONS----END----------------------
        
        IF (PRINT_EACH_PART .EQV. .FALSE.) THEN
            CALL write_2DR(ABsynapse,'ABsynapse')
            
            CALL write_2DR(VSave,'Vsave')
            CALL write_2DR(stdp_save,'stdpc')
            CALL write_2DR(Isyn_save,'IsynS')
            
            CALL write_1DR(VEneu,"VEneu")
            CALL write_1DR(VIneu,"VIneu")
            
            CALL write_1DR(gsumee,"Gsumee")
            CALL write_1DR(gsumei,"Gsumei")
            CALL write_1DR(gsumie,"Gsumie")
        END IF
    END DO
    !!-----------------LOOP OVER ENSEMBLE-----END--------------------------
    CALL SYSTEM_CLOCK(c2)
    systime= (c2-c1)/time_rate
    PRINT*, 'Total spend time is = ', systime, ' seconds.'
END PROGRAM LIFpop
