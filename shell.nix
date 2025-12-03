{ pkgs ? import <nixpkgs> {} }:
let
  python = pkgs.python313;
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
