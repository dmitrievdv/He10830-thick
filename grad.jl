module MagGrad

using Printf

function dvdr(r::Real, θ::Real) :: Real
	return r^(-1.5)*cos(θ)^2/√(4-3*sin(θ)^2)
end 

function dvdt(r::Real, θ::Real) :: Real
	return -2*r^(-0.5)*sin(θ)*cos(θ)/√(4-3*sin(θ)^2)*(3*cos(θ)^2/(4-3*sin(θ)^2)-2)
end

function dudt(r::Real, θ::Real) :: Real
	return -r^(-1.5)/√(4-3*sin(θ)^2)*(cos(θ)^2-sin(θ)^2+3*sin(θ)^2*cos(θ)^2/(4-3*sin(θ)^2))
end

function dudr(r::Real, θ::Real) :: Real
	return 3*r^(-2.5)*cos(θ)*sin(θ)/2/√(4-3*sin(θ)^2)
end

function u(r::Real, θ::Real) :: Real
	return -√(1/r)*cos(θ)^2*2/√(4-3*sin(θ)^2)
end

function v(r::Real, θ::Real) :: Real
	return -√(1/r^3)*cos(θ)*sin(θ)/√(4-3*sin(θ)^2)
end
 
function grad(r::Real, θ::Real, α::Real, β::Real) :: Real
	n = [cos(α), 1/r*sin(α)*cos(β), 1/r/sin(θ)*sin(α)*sin(β)]
	dv = [dvdr(r,θ) dvdt(r,θ); dudr(r,θ) dudt(r,θ)]
	return (dv[1,1]*n[1]^2 + r*sin(θ)^2*n[3]^2*v(r,θ) + r^2*dv[2,1]*n[2]*n[1]+
                    dv[1,2]*n[1]*n[2] + r*v(r,θ)*n[2]^2 + r^2*dv[2,2]*n[2]^2+
                    r^2*sin(θ)*cos(θ)*u(r,θ)*n[3]^2)
end

function grad_integrated_all_directions(r :: Real, θ :: Real; n = 100 :: Integer) :: Real
	αarray = [0:π/n:π;]
	βarray = [0:2*π/n:2*π;]
	s = 0
	for α ∈ αarray, β ∈ βarray
		s = s + abs(grad(r, θ, α, β)) * sin(α)*2*π/n*π/n 
	end
	return s
end

function grad_integrated_over_star(r :: Real, θ :: Real; n = 100 :: Integer) :: Real
	ϕ = 0

	try
		ϕ = asin(1/r)
	catch e
		if isa(e, DomainError)
			ϕ = π/2
		end
	end

	αarray = [π-ϕ:ϕ/n:π;]
	βarray = [0:2*π/n:2*π;]
	s = 0
	for α ∈ αarray, β ∈ βarray
		s = s + abs(grad(r, θ, α, β)) * sin(α)*ϕ/n*2*π/n
	end
	return s
end

function get_grid(rmi :: Real, rmo :: Real; acc = (5,20))
	mr, mθ = acc
	grid = Array{Tuple{Float64,Float64}}(undef, mr, mθ)
	rmags = [rmi:(rmo-rmi)/(mr-1):rmo;]
	i = 1
	for rmag ∈ rmags
		j = 1
		θstart = asin(√(1/rmag))
		θend = π/2
		θspan = θend - θstart
		θstep = θspan/mθ
		# θend = θend - θstep
		θs = [θstart:θstep:θend;]
		print(size(θs), "\n")
		for θ ∈ θs
			if j <= mθ
				grid[i,j] = (rmag, θ)
			end
			j+=1
		end 
		i += 1
	end 
	return grid
end

function get_source_on_grid(grid :: Array{T,2}) :: Array{Float64, 2} where T <: Tuple{Real, Real}
	mr, mθ = size(grid)
	source = zeros(Float64, mr, mθ)
	for i ∈ [1:mr;]
		for j ∈ [1:mθ;]
			rmag, θ = grid[i,j]
			r = rmag*sin(θ)^2
			source[i,j] = grad_integrated_over_star(r,θ)/grad_integrated_all_directions(r,θ)
		end 
	end 
	return source
end

function calc_nu(source :: Real; nl = 1e4, λ = 1.0830e4, Ts = 4e3, gu = 5, gl = 3)
	c = 2.99792458e10
  	h = 6.626176e-27
  	k = 1.380649e-16
	ν = c/(λ*1e-8)
	# print(h*ν, "\n")
	nu = 1/(((exp(h*ν/(k*Ts))-1)/source + 1)*gl/gu/nl)
	return nu
end 

function write_source(source :: Array{R, 2}, grid :: Array{T, 2}) where R<:Real where T<: Tuple{Real, Real}
	mr, mθ = size(source)
	out = open("He10830_thick_popul.dat", "w")
	@printf(out, "# %i %i \n", mr, mθ)
	@printf(out, "#\n")
	for i ∈ [1:mr;]
		for j ∈ [1:mθ;]
			rmag, θ = grid[i,j]
			@printf(out, "%.2f %.7f %.4e %.4e %.4e %.4e \n", rmag, θ, 1e4, source[i,j], 1e4, calc_nu(source[i,j]))
		end
		write(out, "\n")
	end
	close(out)
end


end