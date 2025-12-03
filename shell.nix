
{ pkgs ? import (builtins.fetchTarball {
    url = "https://github.com/NixOS/nixpkgs/archive/25.05.tar.gz";
  }) { config.allowUnfree = true; }
}:

let
  python = pkgs.python312;
  pythonPkgs = python.pkgs;
in

pkgs.mkShell {
  buildInputs = [
    python
    pythonPkgs.numpy
    pythonPkgs.sympy
  ];

  # Optional: make Python the default shell
  shellHook = ''
    echo "numpy and sympy activated."
    unset DRI_PRIME
  '';
}
