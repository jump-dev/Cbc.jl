using BinDeps

@BinDeps.setup

function validate(name,handle)
    try
        # Pre 2.8.12 doesn't have this defined
        p = dlsym(handle, :Cbc_setInitialSolution)
        return p != C_NULL
    catch
        return false
    end
end

@unix_only begin
    libclp = library_dependency("libclp",aliases=["libClp"])
    libcbcsolver = library_dependency("libcbcsolver",aliases=["libCbcSolver"], validate=validate)
end
@windows_only begin
    using WinRPM
    push!(WinRPM.sources, "http://download.opensuse.org/repositories/home:/kelman:/mingw-coinor/openSUSE_13.1")
    libclp = library_dependency("libclp",aliases=["libClp-1"])
    libcbcsolver = library_dependency("libcbcsolver",aliases=["libCbcSolver-3"], validate=validate)
    provides(WinRPM.RPM, "coin-or-Cbc", [libclp,libcbcsolver], os = :Windows)
end

cbcname = "Cbc-2.8.12"

provides(Sources, URI("http://www.coin-or.org/download/source/Cbc/$cbcname.tgz"),
    [libclp,libcbcsolver], os = :Unix)

@osx_only begin
    using Homebrew
    provides( Homebrew.HB, "coinmp", [libclp, libcbcsolver], os = :Darwin )
end

prefix=joinpath(BinDeps.depsdir(libclp),"usr")
patchdir=BinDeps.depsdir(libclp)
srcdir = joinpath(BinDeps.depsdir(libclp),"src",cbcname)

provides(SimpleBuild,
    (@build_steps begin
        GetSources(libclp)
        @build_steps begin
            ChangeDirectory(srcdir)
            `./configure --prefix=$prefix --enable-dependency-linking --without-blas --without-lapack`
            `make install`
        end
    end),[libclp,libcbcsolver], os = :Unix)

@BinDeps.install [:libclp => :libclp, :libcbcsolver => :libcbcsolver]
