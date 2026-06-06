# runtests.jl — MATCRO-Julia regression tests (Fortran gfortran golden refs)
using Test, Printf

@info "MATCRO-Julia Regression Tests"

include("../lib/01_constants.jl")
include("../lib/02_io.jl")
include("../lib/03_hour_interpolation.jl")
include("../lib/04_radiation.jl")
include("../lib/05_1_photosynthesis_C3.jl")
include("../lib/05_2_photosynthesis_C4.jl")
include("../lib/05_photosynthesis.jl")
include("../lib/06_crop.jl")
include("../lib/07_soil.jl")

# ===== Fortran golden refs (gen_ref.exe gfortran -O3) =====
const F_TMP = 289.758766; const F_RSD = 441.7420202; const F_WND = 2.0887081
const F_QPSN = 199.1360685; const F_QPSH = 34.1411229
const F_VSN  = 2.021547e-5; const F_VSH  = 2.021014e-5
const F_DVS  = 0.0543970; const F_SLN  = 0.7525160

@testset "MATCRO vs Fortran (Maize)" verbose=true begin

@testset "interpolate_time" begin
    r = interpolate_time(;
        doy=125, prev_doy=124, next_doy=126, hour=16.5,
        lat=39.75, Δt=3600,
        tmax_prev=299.51755, tmax=292.34836, tmax_next=290.21988,
        tmin_prev=284.44312, tmin=275.42578, tmin_next=277.36337,
        radiation=318.52780, precip=2.2926065e-6,
        humidity=0.008, wind=3.0,
        pressure=100000.0, ozone=0.0,
        int_sinb=2.8480672e4, wind_height=10.0,
    )
    @test r.temperature ≈ F_TMP atol=1e-6
    @test r.radiation   ≈ F_RSD atol=1e-6
    @test r.wind        ≈ F_WND atol=1e-6
end

@testset "calc_radiation (Soybean)" begin
    r = calc_radiation(;
        leaf_nitrogen=1.7683645, kn=0.11,
        shortwave_radiation=519.86699, LAI=3.6562051,
        RLFv=0.105, TLFv=0.07, RLFn=0.58, TLFn=0.25,
        lat=38.25, doy=185, hour=15.5,
        crop_type=SOYBEAN, development_stage=0.5033372,
    )
    @test r.PAR_abs_sunlit_leaf ≈ 138.9702547 rtol=1e-6
    @test r.PAR_abs_shade_leaf  ≈ 22.0396626  rtol=1e-6
    @test r.Vmax_sunlit_leaf    ≈ 6.319203e-5  rtol=1e-6
    @test r.Vmax_shade_leaf     ≈ 5.586533e-5  rtol=1e-6
end

@testset "calc_radiation (Maize)" begin
    r = calc_radiation(;
        leaf_nitrogen=0.58305790, kn=0.30,
        shortwave_radiation=F_RSD, LAI=0.0052646415,
        RLFv=0.105, TLFv=0.07, RLFn=0.58, TLFn=0.25,
        lat=39.75, doy=125, hour=16.5,
        crop_type=MAIZE, development_stage=0.056261203,
    )
    @test r.PAR_abs_sunlit_leaf ≈ F_QPSN atol=1e-4
    @test r.PAR_abs_shade_leaf  ≈ F_QPSH atol=1e-4
    @test r.Vmax_sunlit_leaf    ≈ F_VSN  atol=1e-10
    @test r.Vmax_shade_leaf     ≈ F_VSH  atol=1e-10
end

# ===== PHSYN: trace-based regression (trace = verified Julia Maize run) =====
# DOY=185 HOUR=15.5: C4 day, moderate LAI, meaningful GPP (trace line 394)
@testset "phsyn C4 day (DOY=185 15.5h)" begin
    r = calc_photosynthesis(;
        Qp_sunlit=154.30821, Qp_shade=63.875701,
        Vmax25_sunlit=1.8555916e-5, Vmax25_shade=1.7611940e-5,
        LAI_sunlit=0.63199945, LAI_shade=0.21800503,
        leaf_temperature=308.15127, wind_speed=3.0781774,
        specific_humidity=1.2974163e-2, pressure=9.8608172e4,
        co2_ppm=416.0, water_stress=1.0, crop_height=2.0,
        EFFCON=0.05, atheta=0.8, btheta=0.95,
        m_H2O=4.0, b_H2O=0.04,
        crop_type=MAIZE,
    )
    @test r.gpp ≈ 1.6034169e-5 rtol=1e-6
    @test r.rsp ≈ 7.8869989e-7 rtol=1e-6
    @test r.tsp ≈ 5.8266366e-5 rtol=1e-6
end

# DOY=184 HOUR=19.5: C4 night, LAI=0.86, meaningful TSP (trace line 338)
@testset "phsyn C4 night (DOY=184 19.5h)" begin
    r = calc_photosynthesis(;
        Qp_sunlit=0.0, Qp_shade=0.0,
        Vmax25_sunlit=0.0, Vmax25_shade=1.8468103e-5,
        LAI_sunlit=0.0, LAI_shade=0.86476657,
        leaf_temperature=302.86688, wind_speed=2.9378247,
        specific_humidity=1.1945875e-2, pressure=9.8893281e4,
        co2_ppm=416.0, water_stress=1.0, crop_height=2.0,
        EFFCON=0.05, atheta=0.8, btheta=0.95,
        m_H2O=4.0, b_H2O=0.04,
        crop_type=MAIZE,
    )
    @test r.gpp ≈ 0.0 atol=1e-15
    @test r.rsp ≈ 0.0 atol=1e-15
    @test r.tsp ≈ 1.4875370e-5 rtol=1e-6
end

# C3 Soybean: Soybean debug trace (example/csv run, DOY=185)
# Soybean params: EFFCON=0.425, m_H2O=9.0, b_H2O=0.01, atheta=0.98
@testset "phsyn C3 day Soybean (DOY=185 15.5h)" begin
    r = calc_photosynthesis(;
        Qp_sunlit=138.97026, Qp_shade=22.039662,
        Vmax25_sunlit=6.3179470e-5, Vmax25_shade=5.5854227e-5,
        LAI_sunlit=1.2692534, LAI_shade=2.3869517,
        leaf_temperature=304.09708, wind_speed=1.3315513,
        specific_humidity=1.5671596e-2, pressure=9.8404977e4,
        co2_ppm=438.0, water_stress=1.0, crop_height=1.5,
        EFFCON=0.425, atheta=0.98, btheta=0.95,
        m_H2O=9.0, b_H2O=0.01,
        crop_type=SOYBEAN,
    )
    # Fortran golden (REAL*16): GPP=3.161e-5 RSP=4.642e-6 TSP=1.216e-4
    @test r.gpp ≈ 3.161041e-5 rtol=1e-6
    @test r.rsp ≈ 4.642498e-6 rtol=1e-6
    @test r.tsp ≈ 1.215785e-4 rtol=1e-6
end
@testset "phsyn C3 night Soybean (DOY=185 20.5h)" begin
    r = calc_photosynthesis(;
        Qp_sunlit=0.0, Qp_shade=0.0,
        Vmax25_sunlit=0.0, Vmax25_shade=5.8513140e-5,
        LAI_sunlit=0.0, LAI_shade=3.6351229,
        leaf_temperature=298.73491, wind_speed=1.3315513,
        specific_humidity=1.5671596e-2, pressure=9.8404977e4,
        co2_ppm=438.0, water_stress=1.0, crop_height=1.5,
        EFFCON=0.425, atheta=0.98, btheta=0.95,
        m_H2O=9.0, b_H2O=0.01,
        crop_type=SOYBEAN,
    )
    # Fortran golden (REAL*16): GPP=0 RSP=3.327e-6 TSP=5.678e-6
    @test r.gpp ≈ 0.0 atol=1e-15
    @test r.rsp ≈ 3.327100e-6 rtol=1e-6
    @test r.tsp ≈ 5.678334e-6 rtol=1e-6
end

# ===== CROP =====
@testset "crop_step" begin
    c = CropState(); c.is_planted = 1; c.accumulated_thermal_time = 79.550547
    calc_development_stage!(c, F_TMP, 3600, 1469.0, 0, 8.0, 30.0, 42.0, 60.0)
    @test c.development_stage ≈ F_DVS rtol=1e-6
    c2 = CropState(); c2.development_stage = 0.17104762; c2.leaf_nitrogen = 0.75229684
    calc_specific_leaf_nitrogen!(c2, 130.0, 0.0, 0.52, 1.0, 1.675, 0.825, 0.004, 416.0, MAIZE)
    @test c2.leaf_nitrogen ≈ F_SLN atol=1e-6
    c3 = CropState(); judge_planting!(c3, 120, 12.0, 120, 0)
    @test c3.is_planted == 1; @test !c3.has_emerged
    c4 = CropState(); c4.development_stage = 0.02; judge_emergence!(c4)
    @test c4.has_emerged; @test c4.leaf_biomass ≈ 1.0 atol=1e-10
end

@testset "crop_step (Soybean)" begin
    # CRODVS: Fortran: DVS=0.5033372 aGDH=725.3089
    s = CropState(); s.is_planted = 1
    s.accumulated_thermal_time = 725.0
    calc_development_stage!(s, 304.09708, 3600, 1441.0, 0, 10.0, 27.0, 34.0, 60.0)
    @test s.development_stage ≈ 0.5033372 rtol=1e-6
    @test s.accumulated_thermal_time ≈ 725.3089 rtol=1e-6

    # CALSLN: Fortran: SLN=1.7683645 (NFERT=65, CO2=438, SLNX=0.15/0.4/0.659)
    s2 = CropState(); s2.development_stage = 0.503
    calc_specific_leaf_nitrogen!(s2, 65.0, 0.15, 0.4, 0.659, 3.5, 1.3, 0.004, 438.0, SOYBEAN)
    @test s2.leaf_nitrogen ≈ 1.7683645 rtol=1e-6
end

# ===== SOIL =====
@testset "soil (IRR=1, wet)" begin
    r = calc_soil_water(;
        layer_water=[0.52,0.52,0.52,0.52,0.41],
        transpiration=0.16657306, W2SF=2.2926065e-6,
        depth_root=0.6075, is_irrigated=1, Δt=3600,
        soil_type_i=2, temperature=F_TMP,
        pressure=9.9845492e4, wind_speed=3.3463462,
        specific_humidity=3.9854897e-3,
        crop_height=0.21732839, is_planted=1, crop_type=MAIZE,
    )
    @test r.water_stress ≈ 1.0 rtol=1e-6
end
@testset "soil (IRR=0, dry, low ETC)" begin
    r = calc_soil_water(;
        layer_water=[0.30,0.31,0.32,0.33,0.29],
        transpiration=0.01, W2SF=0.0,
        depth_root=0.5, is_irrigated=0, Δt=3600,
        soil_type_i=2, temperature=F_TMP,
        pressure=9.9845492e4, wind_speed=3.3463462,
        specific_humidity=3.9854897e-3,
        crop_height=0.5, is_planted=1, crop_type=MAIZE,
    )
    @test r.water_stress ≈ 0.6575567 rtol=1e-6
end
@testset "soil (IRR=0, mid, high ETC)" begin
    r = calc_soil_water(;
        layer_water=[0.32,0.33,0.34,0.35,0.31],
        transpiration=0.1, W2SF=0.0,
        depth_root=0.5, is_irrigated=0, Δt=3600,
        soil_type_i=2, temperature=F_TMP,
        pressure=9.9845492e4, wind_speed=3.3463462,
        specific_humidity=3.9854897e-3,
        crop_height=0.5, is_planted=1, crop_type=MAIZE,
    )
    @test r.water_stress ≈ 0.9397901 rtol=1e-6
end

end
