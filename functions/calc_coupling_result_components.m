function results = calc_coupling_result_components(scenario)
    % return structure 'obj' with computed components:
    % plate motion: Rmat, Vpl, Vpl_ss, Vpl_ds
    % slip: strikeSlip, dipSlip, rakeSlip, rakePerpSlip, coupling
    % stress rates: strikeStress, dipStress, rakeStress, rakePerpStress
    % misfit characteristics: chi2, abic
    %
    % Eric Lindsey, June 2019

    % misfit characteristics:
    results.expNumber = scenario.expNumber;
    results.chi2 = ((scenario.predVector - scenario.dataVector)'*inv(scenario.dataCovarianceMatrix)*(scenario.predVector - scenario.dataVector))/length(scenario.dataVector);
    
    % Computation of ABIC is slow: user can compute separately if desired
    %obj.abic = abic_alphabeta(scenario);

    %get the R matrix and plate slip rate
    results.Rmat=scenario.sources{1}.Rmat;
    results.Vpl=scenario.sources{1}.Vpl;
    results.Vpl_ss = results.Rmat(1:end/2,1:end/2) * results.Vpl;
    results.Vpl_ds = results.Rmat(end/2+1:end,1:end/2) * results.Vpl;

    % get the model vector in both strike/dip and rake/rake-perp components
    if isfield(scenario.userParams, 'faultOptions')
        if strcmp(scenario.userParams.faultOptions, 'rakeCoordinates')
            m = results.Rmat * scenario.modelVector;
            results.rakeSlip = - scenario.modelVector(1:end/2);
            results.rakePerpSlip = - scenario.modelVector(end/2+1:end);
        elseif strcmp(scenario.userParams.faultOptions, 'rakeFixed')
            m = results.Rmat(:,1:end/2) * scenario.modelVector;
            results.rakeSlip = - scenario.modelVector;
            results.rakePerpSlip = 0*scenario.modelVector;
        end
    else
        m = scenario.modelVector;
        results.rakeSlip = -results.Rmat(:,1:end/2)' * scenario.modelVector;
        results.rakePerpSlip = -results.Rmat(:,end/2+1:end)' * scenario.modelVector;
    end

    % convert to slip rate magnitude and coupling
    results.strikeSlip = m(1:end/2);
    results.dipSlip = m(end/2+1:end);
    
    % for coupling we compare the slip rate in the rake direction only
    results.coupling = results.rakeSlip ./ results.Vpl;
    
    % compute stress kernels, if needed
    if isempty(scenario.sources{1}.KK)
        [scenario.sources{1}.KK, scenario.sources{1}.Ksn, scenario.sources{1}.Kdn]=unicycle_stress_kernel(scenario.sources{1}.geom, scenario.userParams.slipComponents,scenario.userParams.stressKernelFolder);
    end
    stress_ss_ds = scenario.sources{1}.KK * m;
    results.strikeStress = stress_ss_ds(1:end/2);
    results.dipStress = stress_ss_ds(end/2+1:end);

    stress_r_rp  = results.Rmat' * stress_ss_ds;
    results.rakeStress = stress_r_rp(1:end/2);
    results.rakePerpStress = stress_r_rp(end/2+1:end);

    % moment rate (and "magnitude rate", if useful)
    results.slipMagnitudeMeters = ((results.strikeSlip/1e3).^2 + (results.dipSlip/1e3).^2).^0.5;
    [results.momentDeficitRate, results.magnitudePerYear] = get_moment_and_magnitude(scenario.sources{1}.geom, results.slipMagnitudeMeters);
    
    % rates per 100yrs
    [results.momentDeficitRate100yr, results.magnitudePerYear100yr] = get_moment_and_magnitude(scenario.sources{1}.geom, results.slipMagnitudeMeters*100);
    
end