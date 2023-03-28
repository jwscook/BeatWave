using Plots, FFTW, PyCall, Statistics, Plots

defaults = "./scratch/"

dir = length(ARGS) == 0 ? defaults : ARGS

dirs = filter.(isdir, readdir.("./" .* dir, join=true))[1]

dirs = vcat([d .* ["/both/", "/wave1/", "/wave2/", "/null/"] for d in dirs]...)


function subtractvalue!(x, y = x[1])
  x .-= y
end

NT = 2^10
for dir in dirs
  @show dir
  
  sdf = pyimport("sdf")
  h = sdf.read(dir * lpad(0, 5, "0") * ".sdf")
  Ex0 = h.Electric_Field_Ex_averaged.data
  Bz0 = mean(h.Magnetic_Field_Bz_averaged.data)
  NG = length(Ex0)
  F = zeros(NG, NT);
  fEx(h) = h.Electric_Field_Ex_averaged.data
  fEy(h) = h.Electric_Field_Ey_averaged.data
  fEz(h) = h.Electric_Field_Ez_averaged.data
  fBy(h) = h.Magnetic_Field_By_averaged.data
  fBz(h) = h.Magnetic_Field_Bz_averaged.data

  energydensityelectrons = zeros(NT)
  energydensityions = zeros(NT)
  energydensityfields = zeros(NT)
  L = h.Grid_Grid.data[1][end]
  dx = L / NG

  for i in 0:size(F, 2)-1
    h = sdf.read(dir * lpad(i, 5, "0") * ".sdf")
    try
      energydensityelectrons[i+1] = sum(h.Derived_Average_Particle_Energy_Electron.data) / L
      energydensityions[i+1] = sum(h.Derived_Average_Particle_Energy_Ion.data) / L
      energydensityfields[i+1] = h.Total_Field_Energy_in_Simulation__J_.data / L
    catch e
      @warn "Failed to read particle or field energy from $dir at step $i"
    end
  end
  #subtractvalue!(energydensityelectrons)
  #subtractvalue!(energydensityions)
  #subtractvalue!(energydensityfields)
  energydensitytotal = energydensityelectrons + energydensityions + energydensityfields
  E0 = Bz0^2 / 2 / (4e-7π)
  plot(energydensityelectrons / E0, label="electrons")
  plot!(energydensityions / E0, label="ions")
  plot!(energydensityfields / E0, label="fields")
  plot!(energydensitytotal / E0, label="total")
  xlabel!("Time [arb]")
  ylabel!("Energy density change [<B>^2 / 2 mu0]")
  savefig(dir * "energydensities.png")


  plot((energydensityions .- energydensityions[1]) ./ energydensityions[1])
  savefig(dir * "ionenergydensities.png")

  for (ffield, str) in ((fEx,"Ex"), (fEy, "Ey"), (fEz, "Ez"), (fBy, "By"), (fBz, "Bz"))
    for i in 0:size(F, 2)-1
      h = sdf.read(dir * lpad(i, 5, "0") * ".sdf")
      try
        F[:, i+1] .= ffield(h)
      catch e
        @warn "File number $i has no $str"
      end
    end
    
    filt = sin.(((1:size(F,2)) .- 0.5) ./ size(F,2) .* pi)';

    heatmap(F)
    xlabel!("Time [arb]")
    ylabel!("Position [arb]")
    savefig(dir * str * "_TX.png")

    Ffilt = F .* filt
    nplt = 128 # number modes to plot int time / freq
    nplx = 128 # number modes to plot int space / k
    # Z is now positive omega and k, in correct order. Would plot as k vs omega though
    # transpose to plot omega against k, log and abs for contrast
    Z = log10.(abs.(reverse(fftshift(fft(Ffilt .- mean(Ffilt)))', dims=1)))
    heatmap(-nplx:nplx, 0:nplt, Z[(end÷2):(end÷2 + nplt), (end÷2-nplx+1):(end÷2+nplx+1)])
    xlabel!("Wavenumber [pixel number]")
    ylabel!("Frequency [pixel number]")
    savefig(dir * str * "_WK.png")
  end
end



