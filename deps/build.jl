using BinaryProvider # requires BinaryProvider 0.3.0 or later

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))
products = [
    LibraryProduct(prefix, ["libCbcSolver"], :libcbcsolver),
    LibraryProduct(prefix, ["libCbc"], :libCbc),
]

# Download binaries from hosted location
bin_prefix = "https://github.com/JuliaOpt/CbcBuilder/releases/download/v2.9.9-1-static"

# Listing of files generated by BinaryBuilder:
download_info = Dict(
    Linux(:aarch64, libc=:glibc, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/CbcBuilder.v2.9.9.aarch64-linux-gnu-gcc4.tar.gz", "e0ac415550540b8ce9a4a9845a9f0a39094a5cff20d4383bdd74ea0bb60667b9"),
    Linux(:aarch64, libc=:glibc, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/CbcBuilder.v2.9.9.aarch64-linux-gnu-gcc7.tar.gz", "9b6af44a9f793b4506a306f12401368f3a99affaf6b63e58d1979ce2aba432a0"),
    Linux(:aarch64, libc=:glibc, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/CbcBuilder.v2.9.9.aarch64-linux-gnu-gcc8.tar.gz", "56120d67093ed9121fd19fbb8ca2f5c0d340010f4a2f169d75552bbb9ba29d5a"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/CbcBuilder.v2.9.9.arm-linux-gnueabihf-gcc4.tar.gz", "b2a0848ff62a53e5aafd316180013ccc799c76e8a3c06b5f9eb296b250579bda"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/CbcBuilder.v2.9.9.arm-linux-gnueabihf-gcc7.tar.gz", "93a4d825c7ee88d9d2a06f752ff00ff7537bd52450cf5bf0bd666d64730eaf8b"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/CbcBuilder.v2.9.9.arm-linux-gnueabihf-gcc8.tar.gz", "2ac1751fb3c68ba432e55e1609b41a49b490883938f7007227db8adaabef9ce3"),
    Linux(:i686, libc=:glibc, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/CbcBuilder.v2.9.9.i686-linux-gnu-gcc4.tar.gz", "7a819889d9fc55ce5211fe2beed2d537908e50cb115f3a104e455b3ddad1c520"),
    Linux(:i686, libc=:glibc, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/CbcBuilder.v2.9.9.i686-linux-gnu-gcc7.tar.gz", "c1aa088f903ad4f2b74e9f9b781dce4da763c491856b2dbf359cf96a0a2cd3e1"),
    Linux(:i686, libc=:glibc, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/CbcBuilder.v2.9.9.i686-linux-gnu-gcc8.tar.gz", "e3b436e841ae62505da74bf9ee27dce29e708ad0e3d1589cd7e896bc8597fad1"),
    Windows(:i686, compiler_abi=CompilerABI(:gcc6)) => ("$bin_prefix/CbcBuilder.v2.9.9.i686-w64-mingw32-gcc6.tar.gz", "6465642f5fdd5a4f73a0025a0d86877ffd42a8ac9ca6a2cb14610f67a947c40f"),
    Windows(:i686, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/CbcBuilder.v2.9.9.i686-w64-mingw32-gcc7.tar.gz", "258638aeb2f148850561fbb779193b0763313e96987a279c9abaa0082b4215a1"),
    Windows(:i686, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/CbcBuilder.v2.9.9.i686-w64-mingw32-gcc8.tar.gz", "f0b3809551eb84c67ff0e171a9d6dd799d2672c822dcb1ab647749239edb5327"),
    MacOS(:x86_64, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/CbcBuilder.v2.9.9.x86_64-apple-darwin14-gcc4.tar.gz", "146196c302be996f7b52350113e0648399e5eabce040c0c2ddd1c3f19336b31d"),
    MacOS(:x86_64, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/CbcBuilder.v2.9.9.x86_64-apple-darwin14-gcc7.tar.gz", "8ff4f99f95895cd275e08c5b984d3e42b5b2a17d9f0c1ff3541d81a10cedb2d8"),
    MacOS(:x86_64, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/CbcBuilder.v2.9.9.x86_64-apple-darwin14-gcc8.tar.gz", "ceefa19c327af6373b0e0ba0bd0d68c9391c1ced4db78253df55c531cb48476f"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/CbcBuilder.v2.9.9.x86_64-linux-gnu-gcc4.tar.gz", "ac505a632965ccbe299d46b60e716511df5dc022d4d3467df70df00810fc596d"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/CbcBuilder.v2.9.9.x86_64-linux-gnu-gcc7.tar.gz", "9525b49e5749d1a548dbb7b8d87c54a51db01271a56eb79351e9884a7101c0f0"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/CbcBuilder.v2.9.9.x86_64-linux-gnu-gcc8.tar.gz", "4533ef0c2761339fd252f379bd303bdca7ce94b6f16e0ccfd62782aa8fe65dca"),
    Windows(:x86_64, compiler_abi=CompilerABI(:gcc6)) => ("$bin_prefix/CbcBuilder.v2.9.9.x86_64-w64-mingw32-gcc6.tar.gz", "e62085d05c53e4ffbd4408f6c87a075ec235913c9c8484ef7e2d99f29b3f03ee"),
    Windows(:x86_64, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/CbcBuilder.v2.9.9.x86_64-w64-mingw32-gcc7.tar.gz", "fbd5c104d8f0fb2947d1dc99b09402f11eed9a18980791d2cd489fd3fa08c71f"),
    Windows(:x86_64, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/CbcBuilder.v2.9.9.x86_64-w64-mingw32-gcc8.tar.gz", "a1ad00252a88419d9126afdd46fea1d3a46ebd9aab7d1187657e9189a4972856"),
)

dependencies = [
#     "https://github.com/juan-pablo-vielma/CglBuilder/releases/download/v0.59.10-1/build_CglBuilder.v0.59.10.jl",
#     "https://github.com/JuliaOpt/ClpBuilder/releases/download/v1.16.11-1/build_ClpBuilder.v1.16.11.jl",
#     "https://github.com/juan-pablo-vielma/OsiBuilder/releases/download/v0.107.9-1/build_OsiBuilder.v0.107.9.jl",
#     "https://github.com/juan-pablo-vielma/CoinUtilsBuilder/releases/download/v2.10.14-1/build_CoinUtilsBuilder.v2.10.14.jl",
#     "https://github.com/juan-pablo-vielma/COINMumpsBuilder/releases/download/v1.6.0-1/build_COINMumpsBuilder.v1.6.0.jl",
#     "https://github.com/juan-pablo-vielma/COINMetisBuilder/releases/download/v1.3.5-1/build_COINMetisBuilder.v1.3.5.jl",
#     "https://github.com/juan-pablo-vielma/COINLapackBuilder/releases/download/v1.5.6-1/build_COINLapackBuilder.v1.5.6.jl",
#     "https://github.com/juan-pablo-vielma/COINBLASBuilder/releases/download/v1.4.6-1/build_COINBLASBuilder.v1.4.6.jl",
#     "https://github.com/juan-pablo-vielma/ASLBuilder/releases/download/v3.1.0-1/build_ASLBuilder.v3.1.0.jl"
]
                    
# Install unsatisfied or updated dependencies:
unsatisfied = any(!satisfied(p; verbose=verbose) for p in products)

# To fix gcc4 bug in Windows
this_platform = platform_key_abi()
if typeof(this_platform)==Windows && this_platform.compiler_abi.gcc_version == :gcc4
   this_platform = Windows(arch(this_platform), libc=libc(this_platform), compiler_abi=CompilerABI(:gcc6))
end
dl_info = choose_download(download_info, this_platform)
                    
if dl_info === nothing && unsatisfied
    # If we don't have a compatible .tar.gz to download, complain.
    # Alternatively, you could attempt to install from a separate provider,
    # build from source or something even more ambitious here.
    error("Your platform (\"$(Sys.MACHINE)\", parsed as \"$(triplet(platform_key_abi()))\") is not supported by this package!")
end

# If we have a download, and we are unsatisfied (or the version we're
# trying to install is not itself installed) then load it up!
if unsatisfied || !isinstalled(dl_info...; prefix=prefix)
#     for dependency in reverse(dependencies)          # We do not check for already installed dependencies
#        download(dependency,basename(dependency))
#        evalfile(basename(dependency))
#     end   
    # Download and install binaries
    install(dl_info...; prefix=prefix, force=true, verbose=verbose)
end

# Write out a deps.jl file that will contain mappings for our products
write_deps_file(joinpath(@__DIR__, "deps.jl"), products, verbose=verbose)
