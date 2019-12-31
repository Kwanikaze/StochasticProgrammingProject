import time

from gurobipy import *
import numpy as np

def run_benders_multicut(scenarios, seed, TIME_LIMIT):
    import Project.data_generation as d
    Data = d.data(scenarios, seed)

    num_scenarios = Data.num_scenarios
    num_locations = Data.num_locations
    # Build Sets
    T = ['T'+ str(num) for num in list(range(0,num_scenarios))]
    I = ['I'+ str(num) for num in list(range(0,num_locations))]
    J = ['J'+ str(num) for num in list(range(0,20))]
    K = ['K'+ str(num) for num in list(range(0,10))]
    L = Data.Lset
    M = Data.Mset

    a = Data.a_it  # injury severity (i,t)
    b = Data.b_ijt  # distance injury to evac (i,j,t)
    c = Data.c_ikt  # distance injury to hospital (i,k,t)
    d = Data.d_jkmt  # speed (j,k,m,t)
    e = Data.e_it  # total people injured (i,t)
    f = Data.f_lt  # number patient bed req (l,t)
    ecap = Data.ecap_mt  # m transport capacity (m,t)
    hcap = Data.hcap_lt  # hospital bed type capacity (l,t)
    enod = Data.enod  # evac capacity (j,t)
    hnod = Data.hnod  # hospital capacity (k,t)
    u = Data.u  # max num hospital sites to open
    v = Data.v  # max num evac sites to open
    wx = Data.wx_scen  # proportion ground evac req (t)
    #wx = Data.wx_it
    p = Data.p_t  # probability of each scenario

    def ModifyAndSolveSP(t, time_limit):
        #Remove any for t in T
        #Modify objective coefficients
        expr = LinExpr()
        for i in I:
            for j in J:
                for k in K:
                    for l in L:
                        for m in M:
                            expr.add(x[i,j,k,l,m], a[i,t] * ((b[i,j,t] + c[i,k,t]) / d[j,k,m,t]))

        SP.setObjective(expr, GRB.MINIMIZE)

        # Modify constraint rhs
        for i in I:
            EvacAllConstrs[i].rhs = e[i,t]

        for m in M:
            CapEvacConstrs[m].rhs = ecap[m,t]

        for l in L:
            CapHospitalConstrs[l].rhs = hcap[l,t]

        for j in J:
            CapEvacFlowConstrs[j].rhs = enod * ysol[j]

        for k in K:
            CapHospitalFlowConstrs[k].rhs = hnod * zsol[k]

        #Can't re-use constraint since RHS must be real value/constant but x.sum('*','*','*','*','*') isn't
        SP.remove(SP.getConstrByName("GroundTransportProportion"))
        GroundTransportProportionConstr[0] = SP.addConstr((x.sum('*', '*', '*', '*', 'G') >= wx[t] * x.sum('*','*','*','*','*')),'GroundTransportProportion')

        for l in L:
            PatientTypeConstrs[l].rhs = f[l,t]

        # Solve and get the DUAL solution
        SP.Params.TIME_LIMIT = time_limit
        SP.Params.Threads = 1
        SP.update()
        SP.optimize()
        pi_sol = {}
        gamma_sol = {}
        alpha_sol = {}
        lambda_sol = {}
        delta_sol = {}
        sigma_sol = {}
        psi_sol =  {}

        for i in I:
            pi_sol[i] = EvacAllConstrs[i].Pi

        for m in M:
            gamma_sol[m] = CapEvacConstrs[m].Pi

        for l in L:
            alpha_sol[l] = CapHospitalConstrs[l].Pi

        for j in J:
            lambda_sol[j] = CapEvacFlowConstrs[j].Pi

        for k in K:
            delta_sol[k] = CapHospitalFlowConstrs[k].Pi

        sigma_sol = GroundTransportProportionConstr[0].Pi

        for l in L:
            psi_sol[l] = PatientTypeConstrs[l].Pi

        SPobj = SP.objVal

        # Check whether a violated Benders cut is found
        CutFound_SP = False
        if(nsol[t] < SPobj - CutViolationTolerance): # Found Benders cut is violated at the current master solution
            CutFound_SP = True
        #SP.write("out.lp")
        return SP.Status == GRB.OPTIMAL, SPobj, CutFound_SP, pi_sol, gamma_sol, alpha_sol, lambda_sol, delta_sol, sigma_sol, psi_sol


    # In[3]:


    # Build the master problem
    MP = Model("MP")
    MP.Params.outputFlag = 0  # turn off output
    MP.Params.method = -1      # -1 is automatic, 1 - dual simplex

    # First-stage decision variables
    y = MP.addVars(J, vtype=GRB.BINARY, name='y')  # select evac site
    z = MP.addVars(K, vtype=GRB.BINARY, name='z')  # select hospital site

    # Second-stage eta
    n = MP.addVars(T,obj=p, name='n')
    MP.modelSense = GRB.MINIMIZE

    #First-stage decision variable constraints
    MP.addConstr((y.sum('*') == v), name='CapEvacSites')
    MP.addConstr((z.sum('*') == u), name='CapHospitalSites')

    MP.update()


    # Build the sub-problems
    SP = Model("SP")
    # Removed p[t]
    # Second stage decision variables, will change objective coefficients for each scenario, t
    x = SP.addVars(I, J, K, L, M, obj=1, name='x')

    EvacAllConstrs = {}
    CapEvacConstrs = {}
    CapHospitalConstrs = {}
    CapEvacFlowConstrs = {}
    CapHospitalFlowConstrs = {}
    GroundTransportProportionConstr = {}
    PatientTypeConstrs = {}


    # Setting RHS of these constraints inside Bender's loop later, so for now initialize them using 0
    # RHS changes based on scenario t or based on first stage candidate solutions y and z
    for i in I:
        EvacAllConstrs[i] = SP.addConstr((x.sum(i, '*', '*', '*', '*') == 0),'EvacAll'+str(i))

    for m in M:
        CapEvacConstrs[m] = SP.addConstr((x.sum('*', '*', '*', '*', m) <= 0),'CapEvac'+str(m))

    for l in L:
        CapHospitalConstrs[l] = SP.addConstr((x.sum('*', '*', '*', l, '*') <= 0),'CapHospital'+str(l))

    for j in J:
        CapEvacFlowConstrs[j] = SP.addConstr((x.sum('*', j, '*', '*', '*') <= 0),'CapEvacFlow'+str(j))

    for k in K:
        CapHospitalFlowConstrs[k] = SP.addConstr((x.sum('*', '*', k, '*', '*') <= 0),'CapHospitalFlow'+str(k))


    GroundTransportProportionConstr[0] = SP.addConstr((x.sum('*', '*', '*', '*', 'G')
                                                       >=0),'GroundTransportProportion')

    for l in L:
        PatientTypeConstrs[l] = SP.addConstr((x.sum('*', '*', '*', l, '*') == 0),'PatientType'+str(l))

    SP.modelSense = GRB.MINIMIZE
    SP.update()
    SP.Params.outputFlag = 0 # turn off output

    # Bender's Loop
    CutViolationTolerance = 0.0001
    CutFound = True
    NoIters = 0
    NoCuts = 0
    BestUB = GRB.INFINITY
    time_remaining = TIME_LIMIT
    MPobj = -1
    while (CutFound and time_remaining > 0):
        NoIters += 1
        CutFound = False

        # Solve MP, get new candidate y, z, and eta(n)

        MP.update()
        MP.Params.TIME_LIMIT = time_remaining
        MP.Params.Threads = 1
        start_time = time.time()
        MP.optimize()
        time_remaining = time_remaining - (time.time() - start_time)


        MPobj = MP.objVal  # Get MP solution
        # print('MPobj: %g' % MPobj)

        ysol = {}
        for j in J:
            ysol[j] = y[j].x

        zsol = {}
        for k in K:
            zsol[k] = z[k].x

        nsol = {}
        for t in T:
            nsol[t] = n[t].x

        UB = 0  # First stage costs are always 0
        solved_all_sp = False
        for t in T:
            # Instead of 1 million, intelligently choose enod and hnod for Big M constraints
            # Smaller enod and hnod for each scenario,t  as sum(e_it) over all injury locations, i
            enod = sum(e[i, t] for i in I)
            hnod = enod

            start_time = time.time()
            isOptimal, Qvalue, CutFound_t, pi_sol, gamma_sol, alpha_sol, lambda_sol, delta_sol, sigma_sol, psi_sol = ModifyAndSolveSP(
                t, time_remaining)
            time_remaining = time_remaining - (time.time() - start_time)
            if ~isOptimal:
                solved_all_sp = False
            UB += p[t] * Qvalue
            if (CutFound_t):
                NoCuts += 1
                CutFound = True
                # Dual has 7 variables corresponding to the 7 Primal constraints
                # RHS of Primal constraints become coefficients of Dual Objective
                expr = LinExpr(n[t] - quicksum(e[i, t] * pi_sol[i] for i in I)
                               - quicksum(ecap[m, t] * gamma_sol[m] for m in M)
                               - quicksum(hcap[l, t] * alpha_sol[l] for l in L)
                               - quicksum(enod * y[j] * lambda_sol[j] for j in J)
                               - quicksum(hnod * z[k] * delta_sol[k] for k in K)
                               # - wx[t] * sigma_sol
                               - quicksum(f[l, t] * psi_sol[l] for l in L))
                MP.addConstr(expr >= 0)

        if (UB < BestUB and solved_all_sp):
            BestUB = UB


    return [MPobj, 'Unsure', 'Unsure', TIME_LIMIT - time_remaining]





