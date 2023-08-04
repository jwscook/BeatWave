using Plots, FFTW, PyCall, Statistics, Plots, DSP

defaults = "runs"

const dirarg = length(ARGS) == 0 ? defaults : ARGS

const basedirs = filter.(isdir, readdir.("./" .* dirarg, join=true))
@show basedirs

const dirs = vcat([d .* ["/both/", "/wave1/", "/wave2/", "/null/"] for d in basedirs]...)

@show dirs

function subtractvalue!(x, y = x[1])
  x .-= y
end

const ELECTRON_CHARGE = -1.60217663e-19
const ELECTRON_MASS = 9.1093837e-31
const μ₀ = 4π*1e-7
const kb = 1.380649e-23

function makeplots(dir; skip2dplots=false)
  NT = 2^10

  @show dir
  
  sdf = pyimport("sdf")
  h = sdf.read(dir * lpad(0, 5, "0") * ".sdf")
  Ex0 = h.Electric_Field_Ex_averaged.data
  Bz0 = mean(h.Magnetic_Field_Bz_averaged.data)
  n0 = mean(h.Derived_Charge_Density_averaged_Electron.data) / ELECTRON_CHARGE
  Ωd = abs(Bz0 * ELECTRON_CHARGE) / (1836 * 2 * ELECTRON_MASS)
  Va = abs(Bz0) / sqrt(1836 * 2 * ELECTRON_MASS * n0 * μ₀)
  τcd = 2π / Ωd
  NG = length(Ex0)
  fEx(h) = h.Electric_Field_Ex_averaged.data
  fEy(h) = h.Electric_Field_Ey_averaged.data
  fEz(h) = h.Electric_Field_Ez_averaged.data
  fBy(h) = h.Magnetic_Field_By_averaged.data
  fBz(h) = h.Magnetic_Field_Bz_averaged.data
  fne(h) = h.Derived_Charge_Density_averaged_Electron.data / ELECTRON_CHARGE
  fni(h) = h.Derived_Charge_Density_averaged_Ion.data / abs(ELECTRON_CHARGE)
  fTe(h) = h.Derived_Temperature_averaged_Electron.data * kb / abs(ELECTRON_CHARGE)
  fTi(h) = h.Derived_Temperature_averaged_Ion.data * kb / abs(ELECTRON_CHARGE)
  fearlyEx(h) = h.Electric_Field_Ex.data
  fearlyEy(h) = h.Electric_Field_Ey.data
  fearlyEz(h) = h.Electric_Field_Ez.data
  fearlyBy(h) = h.Magnetic_Field_By.data
  fearlyBz(h) = h.Magnetic_Field_Bz.data

  energydensityelectrons = zeros(NT)
  energydensityions = zeros(NT)
  energydensityfields = zeros(NT)
  times = zeros(NT)

  L = h.Grid_Grid.data[1][end]
  dx = L / NG

  for i in 0:NT-1
    h = sdf.read(dir * lpad(i, 5, "0") * ".sdf")
    times[i+1] = h.Header["time"] / τcd
    try
      energydensityelectrons[i+1] = mean(h.Derived_Average_Particle_Energy_averaged_Electron.data) * n0 * L / L
      energydensityions[i+1] = mean(h.Derived_Average_Particle_Energy_averaged_Ion.data) * n0 * L / L
      energydensityfields[i+1] = h.Total_Field_Energy_in_Simulation__J_.data / L
    catch e
      @warn "Failed to read particle or field energy from $dir at step $i"
      throw(e)
    end
  end
  #subtractvalue!(energydensityelectrons)
  #subtractvalue!(energydensityions)
  #subtractvalue!(energydensityfields)
  E0 = Bz0^2 / 2 / (4e-7π)
  energydensityfields .-= E0
  energydensitytotal = energydensityelectrons + energydensityions + energydensityfields
  plot(energydensityelectrons / E0, label="electrons")
  plot!(energydensityions / E0, label="ions")
  plot!(energydensityfields / E0, label="fields")
  plot!(energydensitytotal / E0, label="total")
  xlabel!("Time [arb]")
  ylabel!("Energy density change [<B>^2 / 2 mu0]")
  savefig(dir * "energydensities.png")

  plot((energydensityions .- energydensityions[1]) ./ energydensityions[1])
  savefig(dir * "ionenergydensities.png")
  if !skip2dplots
      loopargs = ((fEx,"Ex", NT, ""), (fEy, "Ey", NT, ""), (fEz, "Ez", NT, ""), (fBy, "By", NT, ""), (fBz, "Bz", NT, ""),
      (fne, "ne", NT, ""),(fni, "Ni", NT, ""),(fTe, "Te", NT, ""),(fTi, "Ti", NT, ""),
      (fearlyEx,"Ex", NT, "early"), (fearlyEy, "Ey", NT, "early"), (fearlyEz, "Ez", NT, "early"), (fearlyBy, "By", NT, "early"), (fearlyBz, "Bz", NT, "early"))
    for (ffield, str, N, handlestr) in loopargs
      work = handlestr == "early" ? zeros(NG, 512) : zeros(NG, NT)
      ts = zeros(size(work, 2))
      for i in 0:size(work, 2)-1
        h = sdf.read(dir * handlestr * lpad(i, 5, "0") * ".sdf")
        ts[i+1] = h.Header["time"]
        try
          work[:, i+1] .= ffield(h)
        catch e
          @warn "File number $i has no $str"
          throw(e)
        end
      end
      
      filt = DSP.tukey(size(work, 2), 0.5)'

      heatmap(work)
      xlabel!("Time [arb]")
      ylabel!("Position [arb]")
      savefig(dir * str * "_" * handlestr * "_TX.png")

      wfilt = work .* filt
      #wfilt .-= mean(wfilt)
      nplt = 256 # number modes to plot int time / freq
      nplx = 1024 # number modes to plot int space / k
      # Z is now positive omega and k, in correct order. Would plot as k vs omega though
      # transpose to plot omega against k, log and abs for contrast
      Z = log10.(abs.(reverse(fftshift(fft(wfilt))', dims=1)))
      X = (-nplx:nplx) .* 2π / L * Va / Ωd
      Y = (0:nplt) .* 2π / maximum(ts) / Ωd
      heatmap(X, Y, Z[(end÷2):(end÷2 + nplt), (end÷2-nplx+1):(end÷2+nplx+1)])
      xlabel!("Wavenumber [Ωd / Va]")
      ylabel!("Frequency [Ωd]")
      savefig(dir * str * "_" * handlestr * "_WK.png")
    end
  end
  return (energydensityelectrons=energydensityelectrons,
          energydensityions=energydensityions,
          energydensityfields=energydensityfields,
          energydensitytotal=energydensitytotal,
          times=times,
          energynormalisation=E0 * 0.01)
end

dict = Dict()
for dir in dirs
  try
    output = makeplots(dir, skip2dplots=false)
    dict[dir] = output
  catch e
    @warn "Error thrown processing $dir with error $e"
    @show e
  end
end
handles = Dict()
handles["null"] = plot()
handles["both"] = plot()
handles["wave1"] = plot()
handles["wave2"] = plot()

for basedir in basedirs
  he = plot()
  hi = plot()
  hf = plot()
  ht = plot()
  for dirtype in ["null", "both", "wave1", "wave2"]
    dir = basedir * "/" * dirtype * "/"
    nt = dict[dir]

    electrons = nt[:energydensityelectrons]
    ions = nt[:energydensityions]
    fields = nt[:energydensityfields]
    total = nt[:energydensitytotal]

    electrons .-= mean(electrons[2:5])
    ions .-= mean(ions[2:5])
    fields .-= mean(fields[2:5])
    total .-= total[1]

    electrons ./= nt[:energynormalisation]
    ions ./= nt[:energynormalisation]
    fields ./= nt[:energynormalisation]
    total ./= nt[:energynormalisation]

    plot!(he, nt[:times], electrons, label=dirtype)
    plot!(hi, nt[:times], ions, label=dirtype)
    plot!(hf, nt[:times], fields, label=dirtype)
    plot!(ht, nt[:times], total, label=dirtype)
  end
  for (h, title) in ((hi, "ions"), (he, "electrons"), (hf, "fields"), (ht, "total"))
    xlabel!(h, "Time [τcd]")
    ylabel!(h, "Energy density change [<B>^2 / 2 μ₀ / 100]")
    stub = split(basedir, "/")[end]
    savefig(h, "$(stub)_$(title).png")
  end
end

