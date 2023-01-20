using Plots, FFTW, PyCall, Statistics


NT = 2^10
for suffix in ("null", "wave1", "wave2", "both")
  dir = "./scratch/" * suffix * "/"
  
  sdf = pyimport("sdf")
  F = zeros(8192, NT);
  fEx(h) = h.Electric_Field_Ex_averaged.data
  fEy(h) = h.Electric_Field_Ey_averaged.data
  fEz(h) = h.Electric_Field_Ez_averaged.data
  fBy(h) = h.Magnetic_Field_By_averaged.data
  fBz(h) = h.Magnetic_Field_Bz_averaged.data

  for (ffield, str) in ((fEx,"Ex"), (fEy, "Ey"), (fEz, "Ez"), (fBy, "By"), (fBz, "Bz"))
    for i in 0:size(F, 2)-1
      h = sdf.read(dir * lpad(i, 5, "0") * ".sdf")
      F[:, i+1] .= ffield(h)
    end
    
    filt = sin.(((1:size(F,2)) .- 0.5) ./ size(F,2) .* pi)';
    heatmap(F[1:1024, 1:1024])
    savefig(dir * str *"_TX_1024_1024.png")
    heatmap(F[1:512, :])
    savefig(dir * str *"_TX_512_$NT.png")

    heatmap(log10.(abs.(reverse(fft(F)', dims=1))))
    savefig(dir * str * "_WK_whole_unfiltered.png")

    #Z = reverse(fftshift(fft((F .- mean(F)) .* filt)[:, 1:end÷2], 1), dims=1)
    Z = reverse(fftshift(fft((F .- mean(F)))[:, 1:end÷2], 1), dims=1)
    # Z is now positive omega and k, in correct order. Would plot as k vs omega though
    # transpose to plot omega against k, log and abs for contrast
    heatmap(log10.(abs.(Z)')[1:end÷2, (end÷4):(3*end÷4)])
    #Z = reverse(fftshift(fft((F .- mean(F)) .* filt)), dims=1)
    #heatmap(log10.(abs.(Z)'))
    savefig(dir * str * "_WK.png")
  end
end



