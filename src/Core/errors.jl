# SPDX-License-Identifier: MPL-2.0
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>

"""
    CompositionError(msg)

Raised when input cannot be treated as a composition: negative or non-finite
entries, an all-zero sample, mismatched dimensions, or zeros under a policy that
refuses them. A refusal is always an exception and never a silently empty result.
"""
struct CompositionError <: Exception
    msg::String
end

Base.showerror(io::IO, e::CompositionError) = print(io, "CompositionError: ", e.msg)
