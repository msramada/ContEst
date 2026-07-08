# ff_rest.jl — SANDBOX: the remaining three examples from ContEst_five_field_examples.md
#   (5) Epidemiology — surveillance + NPI        → H2/LQG-SDP (linearized SEIR about endemic eq.)
#   (2) Cybergenetics — optogenetic gene control → H2/LQG-SDP (θ_h augments A via reporter dynamics)
#   (1) Neuroscience — closed-loop DBS           → eKF–MPC (nonlinear tanh oscillator network)
# Shared resource budget couples actuation θ_f and sensing θ_h. Reports baseline→joint
# and the joint-vs-sequential co-design gate. Does NOT touch the paper or existing code.
#
# Run:  julia --project benchmarks/ff_rest.jl
# ─────────────────────────────────────────────────────────────────────────────
include("../src/eKF.jl"); include("../src/MPCs.jl"); include("../src/Simulate.jl")
include("../src/LQR.jl"); include("../src/BFGS.jl")
using LinearAlgebra, ForwardDiff, Printf, Optim

zoh(Ac,Bcu,dt)=(mm=size(Bcu,2); Ma=exp([Ac Bcu; zeros(mm,size(Ac,1)+mm)].*dt); (Ma[1:size(Ac,1),1:size(Ac,1)], Ma[1:size(Ac,1),size(Ac,1)+1:end]))

# ── generic H2 co-design gate; A may be a Matrix OR a function θ→A ────────────
function gate_h2(name; A, Bfun, W, Q, R, C, Vfun, θnom, θlb, θub, idxf, idxh, cb, α=1.0)
    A0 = A isa Function ? A(θnom) : A; n=size(A0,1); m=size(Bfun(θnom),2); ny=size(C,1)
    lqr=LQR_θ(n,m); filt=filter_θ(n,ny)
    res(θ)=sum(θ); B_res=sum(θnom)
    Jdes(θ)=cb*(res(θ)-B_res)^2; gdes(θ)=fill(2cb*(res(θ)-B_res),length(θ))
    Amodel(t)= A isa Function ? A(t) : A
    det(θ)=(r=lqr(t->(Amodel(t),Bfun(t),W,Q,R),θ); (r[1],r[2]))
    est(θ)=(r=filt(A, _->C, Vfun, W, Q, θ); (r[1],r[2]))     # filter_θ accepts Af matrix or function
    function evalJ(θ)
        Jd,gd=det(θ); Je,ge=est(θ); (Jd+α*Je+Jdes(θ), gd.+α.*ge.+gdes(θ))
    end
    f,g! = cached_objective(evalJ); rel=verify_gradient(f,g!,θnom)[3]
    _,θj,_=multistart_design(f,g!,θlb,θub;n_starts=5,θ_nom=θnom,seed=20240624,g_tol=1e-7,iterations=150)
    embed(θ0,idx,s)=[k in idx ? s[k-first(idx)+1] : θ0[k] for k in eachindex(θ0)]
    sub(fn,θ0,idx)=(r=optimize(s->fn(embed(θ0,idx,s)),θlb[idx],θub[idx],copy(θ0[idx]),Fminbox(BFGS()),Optim.Options(iterations=200)); embed(θ0,idx,r.minimizer))
    θh=sub(θ->α*est(θ)[1]+Jdes(θ), θnom, idxh)
    θsq=sub(θ->det(θ)[1]+Jdes(θ), θh, idxf)
    parts(θ)=(Jd=det(θ)[1],Je=est(θ)[1],Jt=det(θ)[1]+α*est(θ)[1]+Jdes(θ))
    bp,jp,sp=parts(θnom),parts(θj),parts(θsq)
    println("\n── ($name) H2/LQG-SDP ──  gradient rel $(round(rel,sigdigits=2))")
    @printf("  J_est %.3f→%.3f (%+.1f%%)   J_det %.3f→%.3f (%+.1f%%)   J_tot %.3f→%.3f (%+.1f%%)\n",
            bp.Je,jp.Je,100*(bp.Je-jp.Je)/bp.Je, bp.Jd,jp.Jd,100*(bp.Jd-jp.Jd)/bp.Jd, bp.Jt,jp.Jt,100*(bp.Jt-jp.Jt)/bp.Jt)
    @printf("  CO-DESIGN GATE: sequential %.3f → joint %.3f = %+.2f%%\n", sp.Jt,jp.Jt,100*(sp.Jt-jp.Jt)/sp.Jt)
    @printf("  θ_f %s→%s ; θ_h %s→%s ; %s\n", round.(θnom[idxf],digits=2),round.(θj[idxf],digits=2),
            round.(θnom[idxh],digits=2),round.(θj[idxh],digits=2),
            (any(abs.(θj.-θlb).<1e-3)||any(abs.(θj.-θub).<1e-3)) ? "CORNER" : "INTERIOR")
end

# ── (5) EPIDEMIOLOGY: linearized 2-region SEIR about an endemic operating point ─
function epidemiology()
    β=0.9; σ=0.3; γ=0.25; mob=0.05; s_=0.35; i_=0.04     # endemic-ish operating fractions
    blk=[-β*i_ 0.0 -β*s_; β*i_ -σ β*s_; 0.0 σ -γ]         # per-region [S,E,I] Jacobian
    Ac=zeros(6,6); Ac[1:3,1:3]=blk; Ac[4:6,4:6]=blk
    for c in 0:2                                          # symmetric mobility coupling (diffusion)
        Ac[1+c,1+c]-=mob; Ac[4+c,4+c]-=mob; Ac[1+c,4+c]+=mob; Ac[4+c,1+c]+=mob
    end
    dt=0.5; Bcu=zeros(6,2); Bcu[3,1]=-i_; Bcu[6,2]=-i_    # intervention lowers I in its region
    (Ad,Bdu)=zoh(Ac,Bcu,dt)
    W=Matrix(Diagonal([1e-4,2e-3,4e-3, 1e-4,2e-3,4e-3]))  # outbreak process noise (E,I)
    Q=Matrix(Diagonal([0.0,5.0,20.0, 0.0,5.0,20.0]))      # penalize prevalence I (and exposed E)
    R=Matrix(0.2*I,2,2); C=[0 0 1.0 0 0 0; 0 0 0 0 0 1.0] # test → report I per region
    v0=0.05
    gate_h2("Epidemiology / SEIR"; A=Ad, Bfun=θ->Bdu*Diagonal(θ[1:2]), W=W, Q=Q, R=R, C=C,
        Vfun=θ->Matrix(Diagonal([v0/θ[3], v0/θ[4]])), θnom=ones(4), θlb=fill(0.3,4), θub=fill(4.0,4),
        idxf=1:2, idxh=3:4, cb=4.0)
end

# ── (2) CYBERGENETICS: linearized optogenetic circuit; θ_h enters A (reporter) ─
function cybergenetics()
    δm=0.2; δp=0.1; δr=0.1; kmat=0.3; α_=1.0; β_=1.0; βr_=1.0; φ_=1.0
    # x=[m,p,r1,r2]; θ=[θα,θβ (θ_f); θβr,θmat,θφ (θ_h)]
    Afun(θ)=[ -δm 0.0 0.0 0.0;
              θ[2]*β_ -δp 0.0 0.0;
              θ[3]*βr_ 0.0 -(δr+θ[4]*kmat) 0.0;
              0.0 0.0 θ[4]*kmat -δr]
    Bfun(θ)=reshape([θ[1]*α_,0.0,0.0,0.0],4,1)
    W=Matrix(Diagonal([2e-3,2e-3,1e-4,1e-4])); Q=Matrix(Diagonal([0.0,10.0,0.0,0.0]))
    R=reshape([0.2],1,1); C=[0.0 0.0 0.0 φ_]; v0=0.05
    # forward-Euler discretization (ForwardDiff-clean; gene dynamics are slow, dt·δ≪1)
    dt=0.25
    Ad(θ)=Matrix(1.0I,4,4) .+ dt.*Afun(θ)
    Bd(θ)=dt.*Bfun(θ)
    gate_h2("Cybergenetics"; A=Ad, Bfun=Bd, W=W, Q=Q, R=R, C=C,
        Vfun=θ->reshape([v0/θ[5]],1,1), θnom=ones(5), θlb=fill(0.3,5), θub=fill(4.0,5),
        idxf=1:2, idxh=3:5, cb=3.0)
end

# ── (1) DBS: 2-population damped oscillators, tanh coupling — eKF–MPC ──────────
function dbs()
    NP=2; dt=0.05; Ω=[1.0,1.3]; γo=[0.03,0.03]; η=1.5; κc=0.55; bstim=[1.0,1.0]
    idxv(i)=2i-1; idxvd(i)=2i
    function f_dbs(x,u,θ)
        v=[x[idxv(i)] for i in 1:NP]; vd=[x[idxvd(i)] for i in 1:NP]; e=θ[1:NP]
        vddot(i)= -2γo[i]*vd[i] - Ω[i]^2*v[i] + sum(κc*tanh(η*v[j]) for j in 1:NP if j!=i) + e[i]*bstim[i]*u[i]
        vcat(([v[i]+dt*vd[i], vd[i]+dt*vddot(i)] for i in 1:NP)...)   # interleaved [v_i, vd_i]
    end
    h_dbs(x,::AbstractVector,::AbstractVector)=[x[idxv(i)] for i in 1:NP]     # LFP per contact
    v0=0.2; εp=0.05
    Σw=Matrix(Diagonal(repeat([1e-5,2e-3],NP)))
    Covars(θ)=(Σw, Matrix(Diagonal([v0/(εp+θ[NP+i]^2) for i in 1:NP])))
    Q=Matrix(Diagonal(repeat([10.0,0.0],NP))); R=Matrix(0.5*I,NP,NP)
    detc(x,u)=x'*Q*x+u'*R*u
    stoch(ζ,u)=(xk=ζ[1:2NP]; Σk=vec_to_mat(ζ[2NP+1:end]); xk'*Q*xk+u'*R*u+tr(Q*Σk))
    umax=2.0; con(::AbstractVector,u)=vcat([u[i]-umax for i in 1:NP],[-u[i]-umax for i in 1:NP])
    N=10; prob=ControlModel(f_dbs,h_dbs,Covars,2NP,NP,NP,N,detc,stoch,con)
    mpc=nonlinear_mpc_θ(prob); ce=certainty_equivalence_mpc(prob)
    x_ic=Float64[0.6,0.0,-0.5,0.0]; Σic=Matrix(Diagonal(repeat([4e-2,4e-2],NP))); ulin=zeros(NP)
    θnom=[1.0,1.0,1.0,1.0]; θlb=fill(0.3,4); θub=[fill(3.0,NP);fill(3.0,NP)]
    idxf=1:NP; idxh=NP+1:2NP; B_res=sum(θnom); cb=6.0
    Jdes(θ)=cb*(sum(θ)-B_res)^2; gJ(θ)=fill(2cb*(sum(θ)-B_res),4)
    function jfg(θ)
        _,Jc,g,_ = mpc(x_ic,Σic,ulin,θ;grad=true)
        (Jc+Jdes(θ), g.+gJ(θ))
    end
    function est_cost(θ)
        W,V=resolve_covars(Covars,θ); F=∇ₓf(zeros(2NP),zeros(NP),θ,f_dbs); H=∇ₓh(zeros(2NP),zeros(NP),θ,h_dbs)
        Σ=Matrix{eltype(θ)}(Σic); J=zero(eltype(θ))
        for _ in 1:N; Σp=F*Σ*F'+W; S=H*Σp*H'+V+1e-6I; K=Σp*H'/S; Σ=(I-K*H)*Σp; J+=tr(Q*Σ); end; J
    end
    efg(θ)=(est_cost(θ)+Jdes(θ), ForwardDiff.gradient(est_cost,θ).+gJ(θ))
    embed(θ0,idx,s)=[k in idx ? s[k-first(idx)+1] : θ0[k] for k in eachindex(θ0)]
    subv(vg,θ0,idx)=(r=optimize(s->vg(embed(θ0,idx,s))[1],(G,s)->(G.=vg(embed(θ0,idx,s))[2][idx]),θlb[idx],θub[idx],copy(θ0[idx]),Fminbox(BFGS()),Optim.Options(iterations=50,g_tol=1e-5)); embed(θ0,idx,r.minimizer))
    jf,jg=cached_objective(jfg); rel=verify_gradient(jf,jg,θnom)[3]
    _,θj,_=multistart_design(jf,jg,θlb,θub;n_starts=3,θ_nom=θnom,seed=20240624,g_tol=1e-5,iterations=60)
    θh1=subv(efg,θnom,idxh); θsq=subv(jfg,θh1,idxf)
    ci(xe,Σe,θ)=mpc(xe,Σe,ulin,θ;grad=false)[1]; cm(xe,Σe,θ)=ce(xe,ulin,θ)
    sc(θ;role)= (r=simulate_mc(prob, role==:b ? cm : ci, θ, Q, x_ic, Σic; n_samples=120, seed=20240624); (Je=r.J_est,Jd=r.J_det,Jt=r.J_est+r.J_det+Jdes(θ)))
    b=sc(θnom;role=:b); j=sc(θj;role=:o); s=sc(θsq;role=:o)
    println("\n── (Neuroscience / closed-loop DBS) eKF–MPC ──  gradient rel $(round(rel,sigdigits=2))")
    @printf("  J_est %.3f→%.3f (%+.1f%%)   J_det %.3f→%.3f (%+.1f%%)   J_tot %.3f→%.3f (%+.1f%%)\n",
            b.Je,j.Je,100*(b.Je-j.Je)/b.Je, b.Jd,j.Jd,100*(b.Jd-j.Jd)/b.Jd, b.Jt,j.Jt,100*(b.Jt-j.Jt)/b.Jt)
    @printf("  CO-DESIGN GATE: sequential %.3f → joint %.3f = %+.2f%%\n", s.Jt,j.Jt,100*(s.Jt-j.Jt)/s.Jt)
    @printf("  θ_f(stim) %s→%s ; θ_h(record) %s→%s ; %s\n", round.(θnom[idxf],digits=2),round.(θj[idxf],digits=2),
            round.(θnom[idxh],digits=2),round.(θj[idxh],digits=2),
            (any(abs.(θj.-θlb).<1e-3)||any(abs.(θj.-θub).<1e-3)) ? "CORNER" : "INTERIOR")
end

if abspath(PROGRAM_FILE)==@__FILE__
    println("="^80); println("  ContEst — remaining field examples (epidemiology, cybergenetics, DBS)"); println("="^80)
    epidemiology()
    cybergenetics()
    dbs()
    println("="^80)
end
