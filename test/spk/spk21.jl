
test_dir = artifact"testdata"
DJ2000 = 2451545

function reference_spk1_velocity(g, dt, refvel, Δ, kqmax, kq)
    fc = zeros(length(g) - 1)
    wc = zeros(length(g) - 2)
    w = zeros(length(g) + 2)

    tp = Δ
    mq2 = kqmax - 2
    ks = kqmax - 1

    fc[1] = 1.0
    for j in 1:mq2
        fc[j+1] = tp / g[j]
        wc[j] = Δ / g[j]
        tp = Δ + g[j]
    end

    for j in 1:kqmax
        w[j] = 1 / j
    end

    jx = 0
    ks1 = ks - 1
    while ks >= 2
        jx += 1
        for j in 1:jx
            w[j+ks] = fc[j+1] * w[j+ks1] - wc[j] * w[j+ks]
        end
        ks = ks1
        ks1 -= 1
    end

    for j in 1:jx
        w[j+ks] = fc[j+1] * w[j] - wc[j] * w[j+ks]
    end

    return @inbounds [
        refvel[1] + Δ * sum(dt[j, 1] * w[j] for j in kq[1]:-1:1),
        refvel[2] + Δ * sum(dt[j, 2] * w[j] for j in kq[2]:-1:1),
        refvel[3] + Δ * sum(dt[j, 3] * w[j] for j in kq[3]:-1:1),
    ]
end

@testset "SPK Type 21" verbose=true begin 
    
    # These kernels are tested against SPICE because CALCEPH performs erroneous 
    # computations on SPK types 21 (they differ from those of SPICE)

    # The first kernel has no epoch directories, whereas the second one does
    kernels = [joinpath(test_dir, "spk21_ex1.bsp"), 
               joinpath(test_dir, "spk21_ex2.bsp")]

    yc1 = zeros(3)

    for kernel in kernels

        ephj = EphemerisProvider(kernel);
        furnsh(kernel)

        desc = ephj.files[1].desc[1]
        head = ephj.files[1].seglist.spk1[1].head

        # Center and target bodies 
        cid = Int(desc.cid)
        tid = Int(desc.tid)

        t1j, t2j = desc.tstart, desc.tend

        ep = t1j:1:t2j
        for j in 1:3000
            
            if iseven(j) 
                cid, tid = tid, cid 
            end

            if j == 1 
                # Test initial time
                tj = t1j 
            elseif j == 2
                # Test final time
                tj = t2j 
            elseif j < 100
                # Test values at the directory epochs
                tj = rand(head.epochs)
            elseif j < 500
                # Test directory handling close to the borders
                tj = min(t2j, max(rand(head.epochs) + randn(), t1j))
            else 
                tj = rand(ep)
            end
            tc = tj/86400

            # Test with Julia
            yj1 = ephem_vector3(ephj, cid, tid, tj);
            yj2 = ephem_vector6(ephj, cid, tid, tj);

            # Test with SPICE
            ys1 = spkpos("$tid", tj, "J2000", "NONE", "$cid")[1]
            ys2 = spkezr("$tid", tj, "J2000", "NONE", "$cid")[1]

            @test yj1 ≈ ys1 atol=1e-13 rtol=1e-14
            @test yj2 ≈ ys2 atol=1e-13 rtol=1e-14

            # Test if AUTODIFF works 
            @test D¹(t->ephem_vector3(ephj, cid, tid, t), tj) ≈ yj2[4:end] atol=1e-9 rtol=1e-13

            # TODO: implement the acceleration and jerk. This functions below cannot be tested!
            # D²(t->ephem_vector3(ephj, cid, tid, t), tj)
            # D³(t->ephem_vector3(ephj, cid, tid, t), tj)

            # D¹(t->ephem_vector6(ephj, cid, tid, t), tj)
            # D²(t->ephem_vector6(ephj, cid, tid, t), tj)
            # D³(t->ephem_vector6(ephj, cid, tid, t), tj)

        end

        # Thread-safe testing 
        tj = shuffle(collect(LinRange(t1j, t2j, 200)))

        pos = zeros(3, length(tj))
        for j in eachindex(tj)
            pos[:, j] =  ephem_vector3(ephj, cid, tid, tj[j])
        end
    
        pos_m = zeros(3, length(tj))
        Threads.@threads for j in eachindex(tj)
            pos_m[:, j] = ephem_vector3(ephj, cid, tid, tj[j])
        end

        @test pos ≈ pos_m atol=1e-14 rtol=1e-14

        kclear()

    end

end

@testset "SPK Type 21 velocity coefficient regression" begin
    head = Ephemerides.SPKSegmentHeader1(1, 0, [0.0], 1, 1, 1, 20)
    cache = Ephemerides.SPKSegmentCache1(head)

    cache.g .= 2.0:21.0
    cache.refvel .= (1.0, -2.0, 3.0)
    cache.dt .= reshape(collect(1.0:60.0), 20, 3)
    cache.kqmax = 19
    cache.kq .= (19, 18, 17)

    Δ = 0.25
    Ephemerides.compute_mda_pos_coefficients!(cache, Δ)
    @test_nowarn Ephemerides.compute_mda_vel_coefficients!(cache, Δ)

    expected = reference_spk1_velocity(cache.g, cache.dt, cache.refvel, Δ, cache.kqmax, cache.kq)
    @test Ephemerides.compute_mda_velocity(cache, Δ) ≈ expected atol=1e-14 rtol=1e-14
end
