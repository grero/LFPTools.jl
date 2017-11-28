using Base.Test
using LFPTools
using PyPlot


"""
Create a signal with embedded line noise
"""
function generate_data()
    fs = 30_000.0
    #create a 10 second long signal
    x = 0:1/fs:10.0
    Y = zeros(length(x))
    #random walk
    for i in 2:length(x)
        Y[i] = Y[i-1] + 0.1*randn()
    end
    Y .+= 1.0*sin.(2*pi*50.0*x)
    
    #create a gamma modulation lasting 100ms at random points
    tt = 0.0:1/fs:0.1
    V = 5.0*exp.(-(0.05-tt).^2/(2*0.02^2)).*sin.(2*pi*40.0*tt)
    t0 = Int64[]
    _t0 = rand(1:2*length(V)) 
    for i in 1:length(x)-length(V)
        if i - _t0 >= 2*length(V)
            _t0 = rand(i:i+2*length(V)) 
            push!(t0, _t0)
            Y[_t0:_t0+length(V)-1] .+= V
        end
    end
    x,Y, t0
end

function test()
    srand(1234)
    x, Y, t0 = generate_data()
    Y2,pp = LFPTools.remove_linenoise(Y)
    loss = norm(Y-Y2)
    @test loss ≈ 370.54919665078774
    @test pp[1] ≈ 0.0002496514539056982
    @test pp[2] ≈ 0.956753843020906
    @test pp[3] ≈ -0.03364934764605101
end

@testset "Remove line noise" begin
    test()
end
