function PE = permutation_entropy(signal, m, delay)
    % Calculates permutation entropy
    % signal: input time series
    % m: embedding dimension (typical values: 3-5)
    % delay: time delay (typically 1)
    
    N = length(signal);
    
    % For short signals, use smaller embedding dimension
    if N <= 30 && m > 3
        m = 3; % Adjust embedding dimension for short signals
    end
    
    % Calculate number of possible patterns
    n_patterns = factorial(m);
    
    % Initialize pattern count
    patterns = zeros(1, n_patterns);
    
    % Extract and count patterns
    for i = 1:N-(m-1)*delay
        % Extract values
        pattern = signal(i:delay:i+(m-1)*delay);
        
        % Get ordinal pattern (returns indices of sorted array)
        [~, sorted_indices] = sort(pattern);
        
        % Convert ordinal pattern to single number (1 to n_patterns)
        pattern_num = perm2num(sorted_indices);
        
        % Increment count
        patterns(pattern_num) = patterns(pattern_num) + 1;
    end
    
    % Calculate probabilities
    patterns = patterns / sum(patterns);
    
    % Remove zero probabilities for entropy calculation
    patterns = patterns(patterns > 0);
    
    % Calculate entropy
    PE = -sum(patterns .* log2(patterns));
    
    % Normalize to [0,1]
    PE = PE / log2(n_patterns);
end

function num = perm2num(perm)
    % Converts permutation to single number (1 to factorial(length(perm)))
    m = length(perm);
    num = 1;
    for i = 1:m
        num = num + (perm(i) - 1) * factorial(m - i);
    end
end