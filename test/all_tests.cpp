/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
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
