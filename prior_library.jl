bounds_library = Dict()

# Turbulent kinetic energy parameters
bounds_library[:CᵂwΔ]  = ( 2.0,  20.0)
bounds_library[:Cᵂu★]  = ( 2.0,  20.0)
bounds_library[:Cᴰ⁻]   = ( 0.0,   2.0)
bounds_library[:Cᴰʳ]   = (-1.0,  10.0)
bounds_library[:CᴰRiᶜ] = (-4.0,   4.0)
bounds_library[:CᴰRiʷ] = ( 0.0,   4.0)

# Mixing length parameters
#
#   Recall σ = σ⁻ (1 + σʳ * step(x, c, w))
#
bounds_library[:Cᴷu⁻]  = ( 0.0,  10.0)
bounds_library[:Cᴷc⁻]  = ( 0.0,  10.0)
bounds_library[:Cᴷe⁻]  = ( 0.0,  10.0)
bounds_library[:Cᴷuʳ]  = (-1.0,   2.0)
bounds_library[:Cᴷcʳ]  = (-1.0,   2.0)
bounds_library[:Cᴷeʳ]  = (-1.0,   2.0)
bounds_library[:CᴷRiᶜ] = (-4.0,   4.0)
bounds_library[:CᴷRiʷ] = ( 0.0,   4.0)
bounds_library[:Cᵇu]   = ( 0.0,   4.0)
bounds_library[:Cᵇc]   = ( 0.0,   4.0)
bounds_library[:Cᵇe]   = ( 0.0,   4.0)
bounds_library[:Cˢu]   = ( 0.0,   4.0)
bounds_library[:Cˢc]   = ( 0.0,   4.0)
bounds_library[:Cˢe]   = ( 0.0,   4.0)

bounds_library[:Cᴬˢu]  = ( 0.0,  2.0)
bounds_library[:Cᴬˢc]  = ( 0.0,  2.0)
bounds_library[:Cᴬˢe]  = ( 0.0,  2.0)
bounds_library[:Cᴬu]   = ( 0.0,  0.01)
bounds_library[:Cᴬc]   = (1e-3,  0.1)
bounds_library[:Cᴬe]   = ( 0.0,  0.01)

# Extras
bounds_library[:Cᵇ]    = ( 0.0,   4.0)
bounds_library[:Cˢ]    = ( 0.0,   4.0)
bounds_library[:Cᴰ]    = ( 0.0,   2.0)
bounds_library[:Cᵟu]   = ( 0.0,  10.0)
bounds_library[:Cᵟc]   = ( 0.0,  10.0)
bounds_library[:Cᵟe]   = ( 0.0,  10.0)

prior_library = Dict()

for p in keys(bounds_library)
    prior_library[p] = ScaledLogitNormal(; bounds=bounds_library[p])
end


