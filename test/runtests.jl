using Test
using Fomo

@testset "Fomo.jl" begin
    
    @testset "Backend" begin
        # Test CPU backend creation
        cpu = backend(:cpu)
        @test cpu isa CPUBackend
        
        # Test CUDA check
        if is_cuda_available()
            @info "CUDA available"
        else
            @info "CUDA not available, skipping GPU tests"
        end
    end
    
    @testset "FD Coefficients" begin
        for order in [2, 4, 6, 8]
            coeffs = get_fd_coefficients(order)
            @test length(coeffs) == order ÷ 2
            @test coeffs[1] > 0
        end
    end
    
    @testset "Ricker Wavelet" begin
        f0 = 25.0f0
        dt = 0.001f0
        nt = 1000
        
        wavelet = ricker_wavelet(f0, dt, nt)
        
        @test length(wavelet) == nt
        @test eltype(wavelet) == Float32
        @test maximum(abs.(wavelet)) > 0
    end
    
    @testset "VelocityModel" begin
        nx, nz = 100, 50
        vp = fill(3000.0f0, nx, nz)
        vs = fill(1800.0f0, nx, nz)
        rho = fill(2200.0f0, nx, nz)
        dx, dz = 10.0f0, 10.0f0
        
        model = VelocityModel(vp, vs, rho, dx, dz; name="test")
        
        @test model.nx == nx
        @test model.nz == nz
        @test model.dx == dx
        @test model.dz == dz
    end
    
    @testset "Medium Initialization" begin
        vp = fill(3000.0f0, 100, 50)
        vs = fill(1800.0f0, 100, 50)
        rho = fill(2200.0f0, 100, 50)
        
        be = backend(:cpu)
        nbc = 10
        fd_order = 4
        
        medium = init_medium(vp, vs, rho, 10.0f0, 10.0f0, nbc, fd_order, be)
        
        @test medium.nx == 100 + 2 * nbc
        @test medium.nz == 50 + 2 * nbc
    end
    
    @testset "Wavefield" begin
        nx, nz = 100, 50
        be = backend(:cpu)
        
        wf = Wavefield(nx, nz, be)
        
        @test size(wf.vx) == (nx, nz)
        @test size(wf.vz) == (nx, nz)
        @test size(wf.txx) == (nx, nz)
        
        wf.vx[50, 25] = 1.0f0
        reset!(wf)
        @test wf.vx[50, 25] == 0.0f0
    end
    
    @testset "SimParams" begin
        dt = 0.001f0
        nt = 1000
        dx = 10.0f0
        dz = 10.0f0
        fd_order = 8
        
        params = SimParams(dt, nt, dx, dz, fd_order)
        
        @test params.dt == dt
        @test params.nt == nt
    end
    
    @testset "Simple Simulation (CPU)" begin
        nx, nz = 50, 30
        vp = fill(3000.0f0, nx, nz)
        vs = fill(1800.0f0, nx, nz)
        rho = fill(2200.0f0, nx, nz)
        dx, dz = 20.0f0, 20.0f0
        
        be = backend(:cpu)
        nbc = 10
        fd_order = 4
        
        medium = init_medium(vp, vs, rho, dx, dz, nbc, fd_order, be)
        
        dt = 0.001f0
        nt = 50  # Short for testing
        
        habc = init_habc(medium.nx, medium.nz, nbc, dt, dx, dz, 3000.0f0, be)
        params = SimParams(dt, nt, dx, dz, fd_order)
        fd_coeffs = to_device(get_fd_coefficients(fd_order), be)
        
        wavefield = Wavefield(medium.nx, medium.nz, be)
        
        src_x = Float32[medium.nx * dx / 2]
        src_z = Float32[50.0]
        wavelet = ricker_wavelet(25.0f0, dt, nt)
        
        rec_x = Float32[100.0, 200.0, 300.0]
        rec_z = Float32[20.0, 20.0, 20.0]
        rec = setup_receivers(rec_x, rec_z, medium; type=:vz)
        
        shot_config = MultiShotConfig(src_x, src_z, wavelet)
        
        results = run_shots!(be, wavefield, medium, habc, fd_coeffs,
                            rec, shot_config, params;
                            show_progress=false)
        
        @test length(results) == 1
        @test size(results[1].gather) == (nt, 3)
    end
    
    @testset "VideoConfig" begin
        config = VideoConfig(fields=[:p, :vx], skip=10)
        @test config.skip == 10
        @test :p in config.fields
    end
    
end # testset
