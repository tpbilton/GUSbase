{
  description = "GUSbase";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem
      (system:
        let
          pkgs = import nixpkgs {
            inherit system;
          };

          GUSbase = with pkgs;
            rPackages.buildRPackage {
              name = "GUSbase";
              src = ./.;
              propagatedBuildInputs = with rPackages;
                [ R6 Rdpack ggplot2 doParallel foreach data_table viridis network devtools ];
            };

          R-with-GUSbase = with pkgs;
            rWrapper.override {
              packages = with rPackages;
                [ GUSbase ];
            };
        in
          with pkgs;
          {
            devShells.default = mkShell {
              buildInputs = [ R-with-GUSbase ];
              shellHook = ''
            mkdir -p "$(pwd)/_libs"
            export R_LIBS_USER="$(pwd)/_libs"
          '';
            };

            packages.default = GUSbase;
          });
}
