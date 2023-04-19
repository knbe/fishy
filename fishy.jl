using DelimitedFiles
using StatsPlots
using LaTeXStrings
using Statistics

P_fdu = readdlm("./data/fiducial.dat")
P_Ωmpε = readdlm("./data/omega_m_plus_epsilon.dat")
P_ΩΛpε = readdlm("./data/omega_lambda_plus_epsilon.dat")

const global ε = 0.0001
const global N = size(P_fdu)[1]

const global Ωm = 0.3
const global ΩΛ = 0.7

#function initialize_plot()
#	plot()
#	plot!(
#		xscale = :log10, 
#		yscale = :log10,
#		legend = :false,
#		color = :dimgrey,
#		minorgrid = true,
#		minorgridalpha = .3,
#		size = (1000,1000),
#		titlefontsize=16, 
#		guidefontsize=16,
#		tickfontsize=12,
#		legendfontsize = 16,
#		linewidth = 3,
#		framestyle = :box
#	)
#end

function powerspectrum()
	plot()
	plot!(
		xscale = :log10, 
		yscale = :log10,
		legend = :false,
		color = :dimgrey,
		minorgrid = true,
		minorgridalpha = .3,
		size = (1000,1000),
		titlefontsize=16, 
		guidefontsize=16,
		tickfontsize=12,
		legendfontsize = 16,
		linewidth = 3,
		framestyle = :box
	)

	plot!(
		P_fdu[:,1], 
		P_fdu[:,2],
		xlabel = "log₁₀k (h Mpc⁻¹)", 
		ylabel = "P(k) (h3 Mpc3)",
		legend = :true,
		label = L"$Ω_m=0.3, Ω_Λ=0.7$",
		title = "MATTER POWER SPECTRUM",
	)
#	plot!(
#		P_Ωmpε[:,1], 
#		P_Ωmpε[:,2],
#		label = L"$Ω_m=0.3+ε, Ω_Λ=0.7$",
#	)
#	plot!(
#		P_ΩΛpε[:,1], 
#		P_ΩΛpε[:,2],
#		label = L"$Ω_m, Ω_Λ+ε=0.7$",
#	)

	savefig("./media/43-fiducialmps.png")

	plot!()
end

function Var()
	return var(P_fdu[:,2])
end

function SSR()
	sum = 0.0
	for i in 1:N
		sum += (P_fdu[i,2] - mean(P_fdu[:,2]))^2
	end
	return sum
end

function dP_dΩm()
	sum = 0.0
	for i in 1:N
		sum += (P_Ωmpε[i,2] - P_fdu[i,2]) / ε
	end
	return sum
end

function dP_dΩΛ()
	sum = 0.0
	for i in 1:N
		sum += (P_ΩΛpε[i,2] - P_fdu[i,2]) / ε
	end
	return sum
end

function fisher()
	F_mm = (1 / SSR()) * dP_dΩm() * dP_dΩm()
	F_mΛ = (1 / SSR()) * dP_dΩm() * dP_dΩΛ()
	F_Λm = (1 / SSR()) * dP_dΩΛ() * dP_dΩm()
	F_ΛΛ = (1 / SSR()) * dP_dΩΛ() * dP_dΩΛ()

	F =   [ F_mm F_mΛ ;
		F_Λm F_ΛΛ ]

	return F
end

function fisher2()
	F_mm = (1 / Var()) * dP_dΩm() * dP_dΩm()
	F_mΛ = (1 / Var()) * dP_dΩm() * dP_dΩΛ()
	F_Λm = (1 / Var()) * dP_dΩΛ() * dP_dΩm()
	F_ΛΛ = (1 / Var()) * dP_dΩΛ() * dP_dΩΛ()

	F =   [ F_mm F_mΛ ;
		F_Λm F_ΛΛ ]

	return F
end

function fisher3()
	F_mm = 0.0
	for i in 1:N
		ssr = (P_fdu[i,2] - mean(P_fdu[:,2]))^2
		dpdm = (P_Ωmpε[i,2] - P_fdu[i,2]) / ε
		dpdΛ = (P_ΩΛpε[i,2] - P_fdu[i,2]) / ε
		F_mm += (1/ssr) * dpdm * dpdm
		#println(F_mm)
	end

	F_mΛ = 0.0
	for i in 1:N
		ssr = (P_fdu[i,2] - mean(P_fdu[:,2]))^2
		dpdm = (P_Ωmpε[i,2] - P_fdu[i,2]) / ε
		dpdΛ = (P_ΩΛpε[i,2] - P_fdu[i,2]) / ε
		F_mΛ += (1/ssr) * dpdm * dpdΛ
	end

	F_Λm = 0.0
	for i in 1:N
		ssr = (P_fdu[i,2] - mean(P_fdu[:,2]))^2
		dpdm = (P_Ωmpε[i,2] - P_fdu[i,2]) / ε
		dpdΛ = (P_ΩΛpε[i,2] - P_fdu[i,2]) / ε
		F_Λm += (1/ssr) * dpdΛ * dpdm
	end

	F_ΛΛ = 0.0
	for i in 1:N
		ssr = (P_fdu[i,2] - mean(P_fdu[:,2]))^2
		dpdm = (P_Ωmpε[i,2] - P_fdu[i,2]) / ε
		dpdΛ = (P_ΩΛpε[i,2] - P_fdu[i,2]) / ε
		F_ΛΛ += (1/ssr) * dpdΛ * dpdΛ
	end

	F =   [ F_mm F_mΛ ;
		F_Λm F_ΛΛ ]

	return F
end

function plotfisher()
	plot()
	plot!(
		size = (1000,1000),
		titlefontsize=16, 
		guidefontsize=16,
		tickfontsize=12,
		legendfontsize = 16,
		linewidth = 3,
		framestyle = :box,
		xlabel = L"$Ω_m$",
		ylabel = L"$Ω_Λ$",
		margin = 10Plots.mm,
	)
	covellipse!(
		[0.3, 0.7],
		fisher3() ./ (5.4e4)^3,
		n_std=1,	
		aspect_ratio=1,
		label="1σ",
	)
	covellipse!(
		[0.3, 0.7],
		fisher3() ./ (5.4e4)^3,
		n_std=2,	
		aspect_ratio=1,
		label="2σ",
	)
	plot!(
		[0.3, 0.3], 
		[0.6985, 0.7015], 
		linestyle = :dot,
		linewidth = 2,
		label = L"$Ω_m = 0.3$",
		widen = false,
	)
	plot!(
		[0.2985, 0.3015], 
		[0.7, 0.7], 
		linestyle = :dot,
		linewidth = 2,
		label = L"$Ω_Λ = 0.7$",
		widen = false,
	)

	savefig("./media/fishercontours.png")
	plot!()

		
end
