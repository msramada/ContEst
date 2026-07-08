# ff_linear.jl — SANDBOX: the two LINEAR examples from ContEst_five_field_examples.md
#   (3) Structural — seismic building         → H2 / LQG-SDP  (LQR + dual Kalman Riccati)
#   (4) Fluid — active flow control / drag     → H∞ output-feedback game Riccati
# Each: actuator authority θ_f (in B) vs sensor precision θ_h (in V) under a shared
# budget; report baseline→joint reductions and the joint-vs-sequential co-design gate.
# Self-contained use of src/ solvers; does NOT modify the paper or existing examples.
#
# Run:  julia --project benchmarks/ff_linear.jl
# ─────────────────────────────────────────────────────────────────────────────
include("../src/LQR.jl")      # dare, dlyap, LQR_θ, filter_θ
include("../src/BFGS.jl")     # verify_gradient, multistart_design, cached_objective
include("../src/Hinf.jl")     # Hinf_of_θ (output-feedback H∞)
using LinearAlgebra, Printf, Optim

zoh(Ac,Bcu,dt) = (mm=size(Bcu,2); Ma=exp([Ac Bcu; zeros(mm,size(Ac,1)+mm)].*dt);
                  (Ma[1:size(Ac,1),1:size(Ac,1)], Ma[1:size(Ac,1),size(Ac,1)+1:end]))

# generic co-design gate driver for an H2 (LQR+filter) model with shared budget
function gate_h2(name; A, Bfun, W, Q, R, C, Vfun, θnom, θlb, θub, idxf, idxh, cb, α=1.0, extra=nothing)
    n=size(A,1); m=size(Bfun(θnom),2); ny=size(C,1)
    lqr=LQR_θ(n,m); filt=filter_θ(n,ny)
    B_res=sum(θnom); res(θ)=sum(θ)
    Jdes(θ)=cb*(res(θ)-B_res)^2
    gdes(θ)=fill(2cb*(res(θ)-B_res),length(θ))
    det(θ)=(r=lqr(t->(A,Bfun(t),W,Q,R),θ); (r[1],r[2]))
    est(θ)=(r=filt(A,_->C,Vfun,W,Q,θ); (r[1],r[2]))
    function evalJ(θ)
        Jd,gd = det(θ); Je,ge = est(θ)
        (Jd + α*Je + Jdes(θ), gd .+ α.*ge .+ gdes(θ))
    end
    f, g! = cached_objective(evalJ)
    rel=verify_gradient(f,g!,θnom)[3]
    _,θj,_=multistart_design(f,g!,θlb,θub;n_starts=5,θ_nom=θnom,seed=20240624,g_tol=1e-7,iterations=150)
    # sequential: θ_h for estimation, then θ_f for control
    embed(θ0,idx,s)=[k in idx ? s[k-first(idx)+1] : θ0[k] for k in eachindex(θ0)]
    sub(fn,θ0,idx)=(r=optimize(s->fn(embed(θ0,idx,s)),
                    θlb[idx],θub[idx],copy(θ0[idx]),Fminbox(BFGS()),Optim.Options(iterations=200)); embed(θ0,idx,r.minimizer))
    θh=sub(θ->α*est(θ)[1]+Jdes(θ), θnom, idxh)
    θsq=sub(θ->det(θ)[1]+Jdes(θ), θh, idxf)
    parts(θ)=(Jd=det(θ)[1], Je=est(θ)[1], Jt=det(θ)[1]+α*est(θ)[1]+Jdes(θ))
    bp,jp,sp=parts(θnom),parts(θj),parts(θsq)
    println("\n── ($name) H2/LQG-SDP ──   gradient rel $(round(rel,sigdigits=2))")
    @printf("  J_est %.3f→%.3f (%+.1f%%)   J_det %.3f→%.3f (%+.1f%%)   J_tot %.3f→%.3f (%+.1f%%)\n",
            bp.Je,jp.Je,100*(bp.Je-jp.Je)/bp.Je, bp.Jd,jp.Jd,100*(bp.Jd-jp.Jd)/bp.Jd, bp.Jt,jp.Jt,100*(bp.Jt-jp.Jt)/bp.Jt)
    @printf("  CO-DESIGN GATE: sequential %.3f → joint %.3f = %+.2f%%\n", sp.Jt, jp.Jt, 100*(sp.Jt-jp.Jt)/sp.Jt)
    @printf("  θ_f: %s → %s ;  θ_h: %s → %s\n",
            round.(θnom[idxf],digits=2), round.(θj[idxf],digits=2), round.(θnom[idxh],digits=2), round.(θj[idxh],digits=2))
    interior = !(any(abs.(θj.-θlb).<1e-3)||any(abs.(θj.-θub).<1e-3))
    println("  optimum: ", interior ? "INTERIOR" : "CORNER")
    extra===nothing || extra(θnom,θj)
    (bp=bp,jp=jp,sp=sp,θj=θj,interior=interior)
end

# ── (3) STRUCTURAL: 3-storey shear building, per-story braces + 3 floor sensors ─
function structural()
    nf=3; κ=280.0; M=Matrix(1.0I,nf,nf)
    K=κ.*[2 -1 0;-1 2 -1;0 -1 1]; Cd=0.005.*K+0.02.*M
    Ac=[zeros(nf,nf) Matrix(1.0I,nf,nf); -inv(M)*K -inv(M)*Cd]; n=2nf; dt=0.01
    Td=[1.0 0 0;-1 1 0;0 -1 1]; Tdf=[Td zeros(nf,nf)]
    Γ=[1.0 -1 0;0 1 -1;0 0 1]; (Ad,Bdu)=zoh(Ac,[zeros(nf,nf);inv(M)*Γ],dt)
    eag=[0.0,0,0,1,1,1]; W=1.0.*(eag*eag')+1e-8*Matrix(I,n,n)
    Q=[Td'*(1e3.*Matrix(I,nf,nf))*Td zeros(nf,nf); zeros(nf,nf) 1.0*Matrix(I,nf,nf)]
    R=Matrix(1e-2*I,nf,nf); C=[Matrix(1.0I,nf,nf) zeros(nf,nf)]; v0=1e-2
    driftrms(θ)=(P=dare(Ad,θ[1].*Bdu,Q,R); Kk=-(R+(θ[1].*Bdu)'*P*(θ[1].*Bdu))\((θ[1].*Bdu)'*P*Ad);
                 sqrt.(max.(diag(Tdf*dlyap(Ad+(θ[1].*Bdu)*Kk,W)*Tdf'),0.0)))
    gate_h2("Structural / seismic"; A=Ad, Bfun=θ->θ[1].*Bdu, W=W, Q=Q, R=R, C=C,
        Vfun=θ->Matrix(Diagonal([v0/θ[2]^2,v0/θ[3]^2,v0/θ[4]^2])),
        θnom=ones(4), θlb=fill(0.3,4), θub=fill(5.0,4), idxf=1:1, idxh=2:4, cb=6.0,
        extra=(θn,θj)->(@printf("  peak interstory-drift RMS: %.3e → %.3e (%+.1f%%)\n",
                        maximum(driftrms(θn)),maximum(driftrms(θj)),100*(maximum(driftrms(θn))-maximum(driftrms(θj)))/maximum(driftrms(θn)))))
end

# ── (4) FLUID: ROM active flow control, H∞ output feedback (game Riccati) ──────
function flow(; γ²=40.0)
    ω=[1.0,2.5,4.0]; ζ=[0.03,0.02,0.04]                    # lightly-damped ⇒ gust-resonant modes
    Ac=zeros(6,6)
    for i in 1:3; b=2i-1; Ac[b,b+1]=1.0; Ac[b+1,b]=-ω[i]^2; Ac[b+1,b+1]=-2ζ[i]*ω[i]; end
    dt=0.05; Bcu=zeros(6,1); for i in 1:3; Bcu[2i,1]=1.0; end
    (Ad,Bdu)=zoh(Ac,Bcu,dt)
    E=zeros(6); for i in 1:3; E[2i]=1.0; end; W=0.1.*(E*E')+1e-6*Matrix(I,6,6)   # gust input
    Q=Matrix(Diagonal([1.0,0.2,1.0,0.2,1.0,0.2])); R=reshape([0.1],1,1)
    Cbase=[1.0 0 0 0 0 0; 0 0 1.0 0 0 0]                   # wall sensors read modes 1,2 (positions)
    v0=0.02
    of=Hinf_of_θ(6,1,2;γ²=γ²)
    model(θ)=(Ad, θ[1].*Bdu, W, Q, R, Diagonal([θ[2],θ[3]])*Cbase, Matrix(Diagonal([v0/θ[2]^2,v0/θ[3]^2])))
    θnom=ones(3); θlb=fill(0.3,3); θub=fill(5.0,3); B_res=3.0; cb=6.0
    Jdes(θ)=cb*(sum(θ)-B_res)^2; gdes(θ)=fill(2cb*(sum(θ)-B_res),3)
    function fg(θ)
        Jc,g,_,_,_,ok = of(model,θ)
        (Jc + Jdes(θ), g .+ gdes(θ))
    end
    f, g! = cached_objective(fg); rel=verify_gradient(f,g!,θnom)[3]
    _,θj,_=multistart_design(f,g!,θlb,θub;n_starts=5,θ_nom=θnom,seed=20240624,g_tol=1e-6,iterations=120)
    embed(θ0,idx,s)=[k in idx ? s[k-first(idx)+1] : θ0[k] for k in eachindex(θ0)]
    dJ(θ)=(r=of(model,θ); (r[4],r[5]))                     # (J_det, J_est)
    sub(fn,θ0,idx)=(r=optimize(s->fn(embed(θ0,idx,s)),θlb[idx],θub[idx],copy(θ0[idx]),Fminbox(BFGS()),Optim.Options(iterations=120)); embed(θ0,idx,r.minimizer))
    θh=sub(θ->dJ(θ)[2]+Jdes(θ), θnom, 2:3)
    θsq=sub(θ->dJ(θ)[1]+Jdes(θ), θh, 1:1)
    parts(θ)=(r=of(model,θ); (Jd=r[4],Je=r[5],Jt=r[1]+Jdes(θ)))
    bp,jp,sp=parts(θnom),parts(θj),parts(θsq)
    println("\n── (Fluid / flow control) H∞ output-feedback, γ²=$γ² ──  gradient rel $(round(rel,sigdigits=2))")
    @printf("  J_est %.3f→%.3f (%+.1f%%)   J_det %.3f→%.3f (%+.1f%%)   J_tot %.3f→%.3f (%+.1f%%)\n",
            bp.Je,jp.Je,100*(bp.Je-jp.Je)/bp.Je, bp.Jd,jp.Jd,100*(bp.Jd-jp.Jd)/bp.Jd, bp.Jt,jp.Jt,100*(bp.Jt-jp.Jt)/bp.Jt)
    @printf("  CO-DESIGN GATE: sequential %.3f → joint %.3f = %+.2f%%\n", sp.Jt, jp.Jt, 100*(sp.Jt-jp.Jt)/sp.Jt)
    @printf("  θ_f(act) %.2f→%.2f ; θ_h(sensors) %s→%s\n", θnom[1],θj[1], round.(θnom[2:3],digits=2), round.(θj[2:3],digits=2))
    println("  optimum: ", !(any(abs.(θj.-θlb).<1e-3)||any(abs.(θj.-θub).<1e-3)) ? "INTERIOR" : "CORNER")
end

if abspath(PROGRAM_FILE)==@__FILE__
    using ForwardDiff
    println("="^80); println("  ContEst — LINEAR field examples (structural H2-SDP + fluid H∞)"); println("="^80)
    structural()
    flow()
    println("="^80)
end
