#!/usr/bin/env python3.5

import sys
import os
from RepoManagement.BasicUtils import changedDir, getFrameDir


def SetupTimeDelNetSim():

    # get Current Script Directory
    CurrentScriptDir = getFrameDir()

    # Initialize Directories relevant to CMake
    CMakeModulesDir = os.path.join(CurrentScriptDir, 'CMakeModules')
    CurrentSourceDir = CurrentScriptDir
    BuildDir = os.path.join(os.getcwd(), 'TimeDelNetSim_build')

    # Add CMakeModulesDir to sys path to get relevant functions
    sys.path.insert(0, CMakeModulesDir)
    from CMakePyHelper import getDefaultPlatformGen
    from CMakePyHelper import CMakeGenCall, CMakeBuildCall

    DefaultGenerator = getDefaultPlatformGen()

    # generate and build
    if not os.path.isdir(BuildDir):
        os.mkdir(BuildDir)
    with changedDir(BuildDir):
        isCMakeGenSuccess = CMakeGenCall(CurrentSourceDir,
                                         Generator=DefaultGenerator,
                                         BuildConfig='Release',
                                         Silent=True)
    if isCMakeGenSuccess:
        CMakeBuildSuccess = CMakeBuildCall(BuildDir,
                                           Target='install',
                                           BuildConfig='Release')

    return CMakeBuildSuccess


if __name__ == "__main__":
    SetupTimeDelNetSim()
