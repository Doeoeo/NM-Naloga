"""
	T = PasDiag(sp, d, zg)
 
Podatkovni tip pasovno diagonalno dominantne matrike, ki ima neniclne elmente le ob diagonali.
P = PasDiag([
    100 2 0 0 0;
    4 400 4 0 0; 
    5 6 700 8 0; 
    0 6 7 200 1;
    0 0 8 5 100
])
"""
struct PasDiag 	        # Tip vektorja naceloma ga sam prepozna
    data::Dict{Int, Vector{Number}}
    size::Int
    nSize::Int
    pSize::Int

    function PasDiag(m::Matrix{Int})
        PasDiag(Float64.(m))
    end
    function PasDiag(m::Matrix{Float64})
        # Main dictionary structure to hold diagonals
        d = Dict{Int, Vector{Int}}()
        nrow = size(m)[1]
        ncol = size(m)[2]
        step = nrow + 1
        nSize = -1
        pSize = 0
        # Check if the matrix is diagonally dominant
        for i = 1:ncol
            if sum(m[i, :]) - 2 * m[i, i] >= 0
                error("Matrix m is not diagonally dominant")
                return nothing
            end
        end

        # Save bottom diagonals and main diagonal to negative indexes 0 to -diagNo
        for i = 1:nrow
            diag = m[i:step:(nrow - i + 1) * nrow]
            if all(==(0), diag)
                break
            end
            d[-i + 1] = diag
            nSize += 1
        end
        # Save top diagonals in positive indexes
        j = 1
        for i = 1 + ncol:ncol:length(m)    
            diag = m[i:step:length(m)]
            println(diag)
            if all(==(0), diag)
                break
            end
            d[j] = diag
            pSize += 1
            j += 1
        end
        new(d, nrow, nSize, pSize)        
    end

    function PasDiag(d::Dict{Int, Vector{Number}}, s::Int, n::Int, p::Int)
        new(d, s, n, p)
    end
end

# Overrides
import Base.size, Base.length, Base.*, Base.Matrix, Base.getindex, Base.setindex!, Base.\
"""
    size(P)

Vrni dimenzije pasovne matrike matrike P
"""
size(P::PasDiag) = P.size

"""
    P[x, y]

Vrni vrednost matrike P na indeksu (x, y)
"""
function getindex(P::PasDiag, I::Vararg{Int, 2})
#function getindex(P::PasDiag, I::Integer)
    # Check if index is out of bounds
    if I[1] > P.size || I[2] > P.size 
        error("Index out of bounds at ", I)
        return 0
    end

    # Find the corresponding diagonal by comparing x and y indices
    diagonal = get(P.data, I[2] - I[1], 0)
    if diagonal == 0; return 0 end

    # Return the correct diagonal element based on x and y
    return diagonal[I[1] >= I[2] ? I[2] : I[1]]
end


"""
    P[x, y] = v

Nastavi vrednost v v matriko P na indeks (x, y)
"""
#check types!
function setindex!(P::PasDiag, v::Number,  I::Vararg{Int, 2})
    # Check if index is out of bounds
    if I[1] > P.size || I[2] > P.size 
        error("Index out of bounds at ", I)
    end

    # Find the corresponding diagonal by comparing x and y indices
    P.data[I[2] - I[1]][I[1] >= I[2] ? I[2] : I[1]] = v
end


"""
    y = P*x

Izracunaj produkt pasovno diagonalne matrike P z vektorjem x
"""
function *(P::PasDiag, x::Vector{Int})
    # Check for correct vector length
    if P.size != length(x)
        error("Incompatible vector length")
        return nothing
    end
    y = zeros(length(x))

    # Sum every dim of x
    for i = 1:length(x)
        # Go over all pos diagonals
        for j = 1:min(P.pSize, P.size - i)
            y[i] += get(P.data, j, 0)[i] * x[i]
        end 
        # Go over all neg diagonals
        for j = max(-P.nSize, 1 - i):-1
            y[i] += get(P.data, j, 0)[i + j] * x[i]
        end
        y[i] += P[i, i] * x[i]
    end
    return y
end

# Converter to handle an integer vector b
function \(P::PasDiag, b::Vector{Int})
    b = Float64.(b)
    return P \ b
end

"""
    y = P \ b

Izracunaj y z obratnim vstavljanjem
"""
function \(P::PasDiag, b::Vector{Float64})
    println("yes")
    u, l = lu(P)
    x = deepcopy(b)
    # Compute first value
    x[P.size] = b[P.size] / u.data[0][P.size]
    for i = P.size - 1:-1:1
        line = 0
        for j = 1:min(u.pSize, u.size - i)
            line += u.data[j][i] * x[i + j]
        end
        x[i] = (b[i] - line) / u.data[0][i]
    end

    return x
end


"""
    L, U = lu(P::PasDiag)

Izračunaj LU razcep brez pivotiranja za pasovno diagonalno matriko matriko `P`. Faktorja `L` in `U` sta 
ravno tako matriki tipa `PasDiag`. 
"""
function lu(P::PasDiag)
    p = deepcopy(P)
    # Go over each col to perform GE
    for i = 2:p.size
        # For each row in the col and compute the ratio
        for j = max(-p.nSize, -P.size + i - 1):-1

            #println("outloop ", "j ", j, " " , p.data[j][i - 1], " ", p.data[0][i - 1])
            #println(p.data)
            p.data[j][i - 1] /= p.data[0][i - 1] 

            # Update each row with the ratio of GE
            # First adjust the positive diagonals because they have a dif index
            # The index of the neg diagonal is the index of the positive diagonal element
            for k = 0:min(p.pSize, p.size - i)
                # Extra stopping condition to avoid pointles 0 multiplications
                # note j is a negative number
                if p.pSize - k + j < 0
                    break
                end
                #println("Pos Loop", j, "  " , k, " ", -1)
                #println("inloop ", " " , p.data[k][-j + i - 1], " ", p.data[j][i - 1], " ", p.data[k + 1][i - 1])
                p.data[k][-j + i - 1] -= p.data[j][i - 1] * p.data[k + 1][i - 1]
            end

            # Now adjust the possible negative diagonals
            for k = max(j + 1, -p.nSize):-1
                #println("Neg Loop", j, "  " , k, " ", -1)
                
                if k - j > p.pSize
                    break
                end
                #println("inloop ", "k ", k, " " , p.data[k][- j + k + 1], " ", p.data[j][i - 1], " ", p.data[k - j][i - 1])
                p.data[k][- j + k + i - 1] -= p.data[j][i - 1] * p.data[k - j][i - 1]
            end 
        end
        #println("1 ", p.data[1])
        #println("0 ", p.data[0])
        #println("-1 ", p.data[-1])
        #println("-2 ", p.data[-2])
        #println("-------------------------")
    end
    # Go over the diagonal and pos diagonals


    # Magic :)
    u = Dict( collect(0:p.pSize)[i] => getindex.(Ref(p.data), 0:p.pSize)[i] for i=1:p.pSize+1 )
    l = Dict( collect(-p.nSize:-1)[i] => getindex.(Ref(p.data), -p.nSize:-1)[i] for i=1:p.nSize )
    u = PasDiag(u, p.size, 0, p.pSize)
    l = PasDiag(l, p.size, p.nSize, -1)

    return (u, l)

end