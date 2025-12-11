using KomaMRIFiles

R = 1e-2
r = 5/11*R  
L_tissue = 1e-2 
L_blood = 2e-2
Δx = 8e-4

#POSITIONS
x = -R:Δx:R
y = -R:Δx:R 
z = -L_tissue/2:Δx:L_tissue/2

xx = reshape(x, (length(x),1,1)) 
yy = reshape(y, (1,length(y),1)) 
zz = reshape(z, (1,1,length(z))) 

# Grid
x = 1*xx .+ 0*yy .+ 0*zz
y = 0*xx .+ 1*yy .+ 0*zz
z = 0*xx .+ 0*yy .+ 1*zz

#PHANTOM
⚪(R) =  (x.^2 .+ y.^2 .<= R^2) # circle of radius R

# -------------- Tissue phantom -----------------
ts = Bool.(⚪(R) - ⚪(r))

PD = 1.0
T1 = 1000e-3
T2 = 42e-3

tissue = Phantom(
    name="Tissue",
    x=x[ts],
    y=y[ts],
    z=z[ts],
    ρ=PD.*ones(length(x[ts])),
    T1=T1.*ones(length(x[ts])),
    T2=T2.*ones(length(x[ts]))
)

# -------------- Blood phantom -------------------
    #POSITIONS
    x = -R:Δx:R
    y = -R:Δx:R 
    z = -L_blood/2:Δx:L_blood/2

    xx = reshape(x, (length(x),1,1)) 
    yy = reshape(y, (1,length(y),1)) 
    zz = reshape(z, (1,1,length(z))) 

    # Grid
    x = 1*xx .+ 0*yy .+ 0*zz
    y = 0*xx .+ 1*yy .+ 0*zz
    z = 0*xx .+ 0*yy .+ 1*zz

bl  = Bool.(⚪(r))

PD = 1.0
T1 = 1000e-3
T2 = 42e-3

# Displacements
Nt = 500

dx = dy = zeros(length(z[bl]), Nt)
dz =z[bl] .+ cumsum(L_blood/Nt .+ zeros(1,Nt), dims=2)

spin_reset = dz .> L_blood/2
for i in 1:size(spin_reset, 1)
    idx = findfirst(x -> x == 1, spin_reset[i, :])
    if idx !== nothing
        spin_reset[i, :]  .= 0
        spin_reset[i, idx] = 1 # Se pone a 1 en el nodo SIGUIENTE al salto (ya que en la función reset_spin_hace Constant{Next})
    end
end

dz[dz .> L_blood/2] .-= L_blood
dz .-= z[bl]

v = 4e-2
period = L_blood/v

blood = Phantom(
    name="Blood",
    x=x[bl],
    y=y[bl],
    z=z[bl],
    ρ =PD.*ones(length(x[bl])),
    T1=T1.*ones(length(x[bl])),
    T2=T2.*ones(length(x[bl])),
    motion=flowpath(dx, dy, dz, spin_reset, Periodic(period, 1.0-1e-6))    
)

# ------------- tissue + blood phantom -----------------
phantom = tissue + blood
write_phantom(phantom, "artery.phantom")
