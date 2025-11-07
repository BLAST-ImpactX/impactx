/* Copyright 2022-2023 The Regents of the University of California, through Lawrence
 *           Berkeley National Laboratory (subject to receipt of any required
 *           approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * This file is part of ImpactX.
 *
 * Authors: Axel Huebl, Chad Mitchell, Eric  Stern
 * License: BSD-3-Clause-LBNL
 */
#include "ImpactX.H"

#include <string>


namespace impactx
{
    void ImpactX::call_hook (std::string const & name)
    {
        if (m_hook.count(name) > 0u)
        {
            if (m_hook[name] != nullptr)
            {
                BL_PROFILE("impactx::hook::" + name);
                m_hook[name](this);
            }
        }
    }

} // namespace impactx
