
function wrapper(script_name)

    % ---------------------------------------------------------------------
    % Fix GPU/Mesa issues:
    % Unset DRI_PRIME because Mesa errors out when it is set to 0 on Wayland
    % ---------------------------------------------------------------------
    try
        setenv("DRI_PRIME", "");   % Clear variable
    catch
        % ignore if setenv fails (very rare)
    end

    % ---------------------------------------------------------------------
    % Force FLTK graphics (more stable than qt on Wayland)
    % ---------------------------------------------------------------------
    avail = available_graphics_toolkits;
    if any(strcmpi(avail, "fltk"))
        graphics_toolkit("fltk");
    else
        warning("FLTK not available. Using default toolkit.");
    end

    printf("=== MATLAB Compatibility Wrapper for Octave ===\n");

    % ---------------------------------------------------------------------
    % Load symbolic pkg
    % ---------------------------------------------------------------------
    try
        pkg load symbolic;
        fprintf("Loaded symbolic package.\n");
    catch
        error("Symbolic package is required but not found.");
    end

        % ---------------------------------------------------------------------
    % Load control pkg
    % ---------------------------------------------------------------------
    try
        pkg load control;
        fprintf("Loaded control package.\n");
    catch
        error("Control package is required but not found.");
    end

    % ---------------------------------------------------------------------
    % Reduce SymPy / Octave spam
    % ---------------------------------------------------------------------
    warning("off", "Octave:sym:numeric_to_sym");
    warning("off", "Octave:sym:deprecated-function");
    warning("off", "Octave:deprecated-function");
    warning("off", "Octave:legacy-function");
    warning("off", "all");

    % ---------------------------------------------------------------------
    % Provide fallback vpa()
    % ---------------------------------------------------------------------
    if ~exist("vpa", "file")
        function y = vpa(x, varargin)
            y = double(x);
        end
    end

    printf("Running script: %s\n", script_name);

    % ---------------------------------------------------------------------
    % Execute the MATLAB script
    % ---------------------------------------------------------------------
    run(script_name);

    printf("=== Script finished ===\n");
    pause;

end

