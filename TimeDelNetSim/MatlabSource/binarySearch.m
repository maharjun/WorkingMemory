%  BINARY SEARCH ALGORITHM
%  Description : 
%     mex function that performs the binary search algorithm to find "item(s)"
%     (the values to be searched for) among some pre-sorted "data" vector.
%     By default, the algorithm returns the index of the first instance of each "item"
%     (if there are multiple copies found), and returns the index of the closest item
%     if the item(s) are not found (although these behaviors can be changed 
%     with appropriate optional input parameters.)
% 
%   Note : by default, the algorithm does not check whether the input data is sorted (since 
%    that would be an O(N) procedure, which would defeat the purpose of the
%    algorithm.  If the input data is not sorted, the output values will be incorrect.
% 
%  
%  Matlab call syntax:
%     pos = binarySearchMatlab(data, items, [dirIfFound], [dirIfNotFound], [checkIfSorted_flag])
% 
%  Matlab compile command:
%     mex binarySearch.c
% 
%  Input: This function requires (pre-sorted) reference data vector "data", 
%   as well as a second input, "items" to search for. "items" can be any size. 
% 
%  Output : "pos" is the same size as the input "items".
%    
%  * Optional input arguments: 'dirIfFound' and 'dirIfNotFound' specify
%    the function's behavior if the items are not found, or if multiple 
%    items are found: (Supply an empty vector [] to leave as the default.)
%    Note, if you like, you can change the default behavior in each case by
%    modifying the DEFAULT values in the #define section below.
% 
%    dirIfFound  specifies the function's behavior if *multiple* copies of the 
%       value in "items" are found.
%      -1, or 'first' : [default] the position of the *first* occurence of 'item' is returned
%      +1, or 'last'  : the position of the *last* occurence of 'item' is returned. 
%       0, or 'any'   : the position of the first item that the algorithm discovers is found
%            (ie. not necessarily the first or last occurence)
% 
%    dirIfNotFound specifies the behavior if the value in "items" is *not* found.
%         0, or 'exact'   : the value 0 is returned.            
%        -1, or 'down'    : the position of the last item smaller than 'item' is returned.
%        +1, or 'up'      : the position of the first item greater than 'item' is returned.
%         2, or 'closest' : [default], the position of the *closest* item closest to 'item' 
%                           is returned
%         0.5, or 'frac'  : the function returns a *fractional* value, indicating, the 
%                           relative position between the two data items between which 'item' 
%                           would be located if it was in the data vector. 
%                           (eg if you are searching for the number 5 (and "data" starts off 
%                           with [ 2, 3, 4, 7,...], then the algorithm returns 3.333, because 
%                           5 is 1/3 of the way between the 3th and the 4th elements of "data".
% 
%   checkIfSorted_flag
%      By default, this program is set *not* to check that the input data vector is sorted.
%     (although you can change this by setting the defined CHECK_IF_INPUT_SORTED_DEFAULT as 1)
%     However, if you provide a non-empty 5th argument the input data will be checked. 
%     (You might use this, for example, while debugging your  code, and remove it later 
%      to improve performance)
% 
%   Example:
%      data = 1:100;
%      items = [pi, 42, -100]
%      binarySearch(data, items)
%      ans =
%           3    42     1
% 
% 
% 
%   Please send bug reports / comments to :
%   Avi Ziskind
%   avi.ziskind@gmail.com
%   
%   last updated: May 2013.
% 
%   update on 5/2/2013:
%     * fixed a memory leak that occurs if strings are passed as 3rd or 4th arguments (you need to 
%       call mxFree if you use mxArrayToString)
%     * added out-of-bounds check to the binary search core (if item is out of bounds, we can skip 
%       the search altogether)
%     * allow for both single-precision or double-precision inputs.
%       ('data' and 'item' can both be either double or single). This is done by having two copies
%       of the core search function, one that is used if 'data' is double, and one for if 'data' is
%       single (and each 'item' is cast accordingly). If anyone knows of a better way to do this 
%       that isn't too cumbersome and won't add too much overhead, please let me know.


