/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#define CATCH_CONFIG_RUNNER  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch/include/catch.hpp"

bool CHAIN_EXTENSIVE = false;

int main(int argc, char* const argv[])
{
    Catch::Session session;

    int returnCode = session.applyCommandLine(argc, argv, Catch::Session::OnUnusedOptions::Ignore);
    if (returnCode != 0)
        return returnCode;

    for (auto token : session.unusedTokens()) {
        printf("Token: %s\n", token.data.c_str());
        if (token.data == "EXTENSIVE")
            CHAIN_EXTENSIVE = true;
    }

    return session.run();
}
