// decifra_hill.sce — Descriptografia Hill (mod 29) com Gauss-Jordan em Scilab
// Requer: dadosDoProblema.txt no mesmo diretório.

//------------------------
// Utilitários de módulo
//------------------------
function y = pmod(A, m)
    // resto positivo elemento a elemento
    y = modulo(A, m);
    y = y + m*(y<0);
endfunction

function [g, x, y] = egcd(a, b)
    // Algoritmo de Euclides estendido (escalares)
    // Retorna g = mdc(a,b) e x,y tais que ax + by = g
    a = int(a); b = int(b);
    x0=1; y0=0; x1=0; y1=1;
    while b<>0
        q = floor(a/b);
        tmp = a - q*b; a = b; b = tmp;
        tmp = x0 - q*x1; x0 = x1; x1 = tmp;
        tmp = y0 - q*y1; y0 = y1; y1 = tmp;
    end
    g=a; x=x0; y=y0;
endfunction

function inva = inv_mod_scalar(a, m)
    a = pmod(a, m);
    [g, x, _] = egcd(a, m);
    if g<>1 then
        error("Sem inversa: gcd("+string(a)+","+string(m)+") = "+string(g));
    end
    inva = pmod(x, m);
endfunction

//-----------------------------------------
// Gauss–Jordan modular para matriz  K^-1
//-----------------------------------------
function Kinv = inv_mod_matrix(K, m)
    [n, p] = size(K);
    if n<>p then error("K não é quadrada."); end
    // Monta [K | I]
    A = pmod([K , eye(n,n)], m);
    r = 1;
    for c=1:n
        // encontra pivô ≠ 0 em A(r,c)
        piv = 0; piv_row = 0;
        for i=r:n
            if pmod(A(i,c),m)<>0 then piv = A(i,c); piv_row = i; break; end
        end
        if piv_row==0 then
            error("Matriz não-invertível mod "+string(m)+" (coluna "+string(c)+").");
        end
        // troca linha r <-> piv_row
        if piv_row<>r then
            tmp = A(r,:); A(r,:) = A(piv_row,:); A(piv_row,:) = tmp;
        end
        // normaliza linha r (divide pelo pivô modulo m)
        inv_piv = inv_mod_scalar(A(r,c), m);
        A(r,:) = pmod(A(r,:)*inv_piv, m);
        // zera as demais linhas na coluna c
        for i=1:n
            if i<>r then
                fator = A(i,c);
                if fator<>0 then
                    A(i,:) = pmod(A(i,:) - fator*A(r,:), m);
                end
            end
        end
        r = r + 1;
    end
    // agora A = [I | Kinv]
    Kinv = pmod(A(:, n+1:2*n), m);
endfunction

//-----------------------------------------
// Conversão texto <-> números (alfabeto)
//-----------------------------------------
function nums = text_to_nums(s)
    // Normaliza para minúsculas
    s = convstr(s, "l");

    // Remove caracteres de controle comuns
    s = strsubst(s, char(13), ""); // CR (\r)
    s = strsubst(s, char(9),  " "); // TAB -> espaço

    // Remove BOM (U+FEFF ou sequência 239 187 191)
    if length(s) > 0 then
        if ascii(part(s,1)) == 65279 then
            s = part(s, 2:length(s));
        elseif length(s) >= 3 then
            if ascii(part(s,1)) == 239 & ascii(part(s,2)) == 187 & ascii(part(s,3)) == 191 then
                s = part(s, 4:length(s));
            end
        end
    end

    // Alfabeto permitido como STRING
    alphabet = ascii([97:122, 46, 32, 44]); // "abcdefghijklmnopqrstuvwxyz. ,"
    acodes   = ascii(alphabet);             // códigos de cada posição (29 elementos)

    n = length(s);
    nums = zeros(1, n);
    for i = 1:n
        ch  = part(s, i);
        idx = find(acodes == ascii(ch));    // << comparação por código
        if idx == [] then
            error("Caractere inválido '"+ch+"' (código "+string(ascii(ch))+") na posição "+string(i)+".");
        end
        nums(i) = idx(1) - 1; // zero-based
    end
endfunction

function s = nums_to_text(nums)
    alphabet = ascii([97:122, 46, 32, 44]); // "abcdefghijklmnopqrstuvwxyz. ,"
    n = size(nums, "*");
    chars = matrix(" ", 1, n);
    for i=1:n
        v = nums(i);
        if v<0 | v>28 then error("Valor fora do alfabeto: "+string(v)); end
        chars(i) = part(alphabet, v+1);
    end
    s = strcat(chars, "");
endfunction

//-----------------------------------------
// Leitura do arquivo e decifragem
//-----------------------------------------
function main()
    m = 29;

    // Matriz K (10x10) do enunciado
    K = [ 2, 18, 20, 11,  5,  0,  4,  8, 10,  1;
          7, 21,  3, 14, 25, 17, 19, 28,  6, 13;
         10,  4, 16,  9,  2, 22,  1, 27, 12,  5;
         24,  8, 15, 23, 11, 19,  0,  3,  7, 26;
          1, 14, 28,  5, 17,  6, 21, 10,  4, 20;
          9,  0, 11, 22,  7, 13, 25,  2, 16, 18;
         12, 26,  4,  1, 20,  8, 14, 23,  5, 27;
         19,  7, 24, 10,  3, 28, 17,  5, 21,  9;
         22, 13,  6, 16,  0, 27,  8, 11, 15,  2;
          5, 12, 23, 18, 26,  9, 13,  1, 24,  7];

// Resolve o arquivo a partir da pasta do próprio script:
dirscript = get_absolute_file_path("pmod.sci");
if dirscript=="" then dirscript = pwd(); end

// Garante o separador correto e evita duplicação:
sep = filesep();
if part(dirscript, length(dirscript))==sep then
    path = dirscript + "dadosDoProblema.txt";
else
    path = dirscript + sep + "dadosDoProblema.txt";
end

    lines = mgetl(path);              // vetor de linhas
    Ctxt  = strcat(lines, "");        // concatena em uma única string

    // Converte para números e valida blocos de 10
    Cnums = text_to_nums(Ctxt);
    n = length(Cnums);
    if modulo(n, 10)<>0 then
        error("O comprimento do texto cifrado ("+string(n)+") não é múltiplo de 10.");
    end

    // Número de blocos de 10 símbolos
    nblocks = n/10;

    // Colunas são blocos (10 x nblocks)
    Cmat = matrix(Cnums, 10, nblocks);

    // Inversa de K (mod 29)
    Kinv = inv_mod_matrix(K, m);

    // *** ATENÇÃO: cifra foi gerada como C = P * K (vetor-linha) ***
    // Então a decodificação correta é P = C * K^{-1} (à direita).
    // Como Cmat está 10xT com blocos em colunas, transponha, multiplique à direita e destroque:
    Pmat = pmod( (Cmat') * Kinv , m )';

    // Volta para texto
    Pnums = matrix(Pmat, 1, n);
    Ptxt  = nums_to_text(Pnums);

    // Exibe resultado
    disp("===== TEXTO DECIFRADO =====");
    disp(Ptxt);
    disp("===========================");
endfunction

// Executa
main();
