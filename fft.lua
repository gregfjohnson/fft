-----------------------------------------------------------------------------
-- Copyright (c) Greg Johnson, Gnu Public Licence v. 2.0.
-----------------------------------------------------------------------------
-- Version 0.1.0.
--[[
    lua fft implementation.
    fft and inverse fft preserve vector lengths.

    you will also need https://github.com/gregfjohnson/cmath.

    Usage:
    -- require returns a function:
    fft = require 'fft'

    vec = {1,2,3}
    result = fft(vec)

    -- do an inverse fft:
    inverseFft = false
    back   = fft(result, inverseFft)

    -- result and back are equal to within tolerance..
    function eq(x, y) return math.abs(x-y) < 1e15 end

    for i = 1, #vec do
        assert(eq(vec[i], back[i]))
    end

    -- complex input vector:
    -- (fft requires cmath internally, but keeps it local.)
    require 'cmath'
    x = fft{1, i, -1, -i}
--]]

local verbose = false
local dbPrint

local cmath  = require 'cmath_anon'
local primes = require 'primes'

local function pvec(v) for i = 1, #v do print(v[i]) end; print() end

local function recFFT(invec, doForward)
    if verbose then print('invec:') pvec(invec) end

    local n = #invec
    local outvec = {}
    for i = 1, #invec do outvec[i] = 0 end

    if #invec == 1 then
        outvec[1] = invec[1]

    else
        local chunkCount = primes.smallestPrimeFactor(n)
        local chunkSize  = n / chunkCount

        local chunks = {}
        for i = 1, chunkCount do chunks[i] = {} end

        for i = 0, n-1 do
            local indexInChunk = i / chunkCount
            local chunkIndex   = i % chunkCount

            table.insert(chunks[chunkIndex + 1], invec[i + 1])
        end

        local results = {}
        for chunk = 1, chunkCount do
            results[chunk] = recFFT(chunks[chunk], doForward)
        end

        local directionalOmega
        if doForward then directionalOmega = cmath.exp(-2 * cmath.pi * cmath.i / n)
        else              directionalOmega = cmath.exp( 2 * cmath.pi * cmath.i / n)
        end

        for i = 0, n-1 do
            outvec[i+1] = results[1][1 + i % chunkSize]
        end

        if verbose then print('initial outvec:') pvec(outvec) end

        for chunk = 1, chunkCount-1 do
            for i = 0, n-1 do
                local omegaPower = directionalOmega ^ (i * chunk)

                if verbose then
                    dbPrint(chunk, chunkSize, i, directionalOmega, omegaPower, results, outvec)
                end

                outvec[i+1] =   outvec[i+1]
                            + results[chunk + 1][i % chunkSize + 1] * omegaPower

                if verbose then print('  -> ', outvec[i+1]) end
            end
        end
    end

    if verbose then print('outvec:') pvec(outvec) end

    return (outvec)
end

local function FFT(invec, forward)
    if forward == nil then forward = true end

    local result = recFFT(invec, forward)

    local len = math.sqrt(#result)
    for i = 1, #result do result[i] = result[i] / len end
    return result
end

function dbPrint(chunk, chunkSize, i, directionalOmega, omegaPower, results, outvec)
    io.write(chunk, ' ', i, ' ',
             tostring(directionalOmega), ' ',
             tostring(omegaPower), ' || ',
             tostring(results[chunk + 1][i % chunkSize + 1]*omegaPower),
             ' ( ', tostring(results[chunk + 1][i % chunkSize + 1]), ' ^= ',
                    tostring(results[chunk + 1][i % chunkSize + 1]*omegaPower), ' ) ',
             ' + ',
             tostring(outvec[i+1]))
end

return FFT
