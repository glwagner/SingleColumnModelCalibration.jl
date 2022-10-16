using Distributions

bounds_library = Dict()

# Turbulent kinetic energy parameters
bounds_library[:CᵂwΔ]  = ( 0.0,  15.0)
bounds_library[:Cᵂu★]  = ( 0.0,  15.0)
bounds_library[:Cᴰ⁻]   = ( 0.0,  4.0)
bounds_library[:Cᴰ⁺]   = ( 0.0,  4.0)
bounds_library[:CᴰRiᶜ] = ( 0.0,  1.0)
bounds_library[:CᴰRiʷ] = ( 0.0,  1.0)

# Mixing length parameters
bounds_library[:Cᴷu⁻]  = ( 0.0, 4.0)
bounds_library[:Cᴷc⁻]  = ( 0.0, 4.0)
bounds_library[:Cᴷe⁻]  = ( 0.0, 4.0)

bounds_library[:Cᴷu⁺]  = ( 0.0, 4.0)
bounds_library[:Cᴷc⁺]  = ( 0.0, 4.0)
bounds_library[:Cᴷe⁺]  = ( 0.0, 4.0)

bounds_library[:CᴷRiᶜ] = ( 0.0, 1.0)
bounds_library[:CᴷRiʷ] = ( 0.0, 1.0)

bounds_library[:Cᵇu]   = ( 0.0, 4.0)
bounds_library[:Cᵇc]   = ( 0.0, 4.0)
bounds_library[:Cᵇe]   = ( 0.0, 4.0)

bounds_library[:Cˢu]   = ( 0.0, 4.0)
bounds_library[:Cˢc]   = ( 0.0, 4.0)
bounds_library[:Cˢe]   = ( 0.0, 4.0)

bounds_library[:Cᵟu]   = ( 0.0, 4.0)
bounds_library[:Cᵟc]   = ( 0.0, 4.0)
bounds_library[:Cᵟe]   = ( 0.0, 4.0)

bounds_library[:Cᴬu]   = ( 0.0,  10.0)
bounds_library[:Cᴬc]   = ( 0.0, 100.0)
bounds_library[:Cᴬe]   = ( 0.0, 100.0)

bounds_library[:Cʰ]    = ( 0.0,   0.1)
bounds_library[:Cʰˢ]   = ( 0.0,   4.0)

bounds_library[:Cʷ★]   = ( 1.0,  4.0)
bounds_library[:Cʷℓ]   = ( 0.0,  4.0)

prior_library = Dict()

for p in keys(bounds_library)
    prior_library[p] = ScaledLogitNormal(; bounds=bounds_library[p])
end


