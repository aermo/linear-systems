(ns
  ^{:author "Antti Rämö",
    :doc ""}
  linear-systems.core)


(defn diag [A] 
  ; (diag A) Returns a diagonal element of a matrix
  ; (diag v) Returns a diagonal matrix with diagonal units v.
  (if (coll? (first A))
    (map nth A (range))
    (map #(concat (repeat %2 0) [%1] (repeat (- (count A) (inc %2)) 0)) A (range))))

(defn matrix-vector-prod [m v] 
  (map #(apply + (map * % v)) m))

(defn matrix-prod [A B]
  (map (fn [ai] (map (fn [bj] (apply + (map * ai bj))) (apply map list B))) A))

(defn pivot
  ([A] (pivot A (range (count A)) 0))
  ([A k] (pivot A (range (count A)) k))
  ([A p k]
    ; Find the row with the greatest value in the column k.
    (let [i (->> 
              ;(map first A)  ; Take first column.
              (map #(nth % k) A)  ; Take column k.
              (map #(* % %))      ; Square piecewise
              (zipmap (range))    ; Add row index.
              (apply max-key val) ; Take the key-value-pair with greatest value.
              (key))              ; Extract the key.
          r (into [i] (remove (partial = i) (range (count A))))] ; New row order.
      {:A (replace (vec A) r)
       :p (replace (vec p) r)
       :r r})))

(defn backward-substitution-sparse [U c] 
; Solves a linear system 'Ux = c' by calculating the backward substitution.
  (if (map? U)
    ; If U is returned from LU-factorization (returns map: {:L L :U U :p p}).
    (backward-substitution-sparse (U :U) c)
    (if (coll? (first c))
      ; If c is a matrix calculate each column separately.
      (map (partial backward-substitution-sparse U) c)
      ; Take upper triangular matrix without zeros and cons vector c.
      (let [U (map cons c U)]
        (reduce
          #(cons 
            (/ (- (first %2) (apply + (map * %1 (nnext %2)))) (second %2))
            %1)
          '() (reverse U))))))


(defn gaussian-elimination
  ([A b] 
    (gaussian-elimination (map #(conj (vec %1) %2) A b)))
  ([A]
    (letfn [(elimination-step [A]
              (map 
                (fn [ai-jj]
                  (map 
                    #(- %1 (/ (* %2 (first ai-jj)) (first (first A))))
                    (rest ai-jj) (rest (first A)))) 
                (rest A)))
            (sub-ge [A k]
              (when-first [a1 A]
                (cons (concat (repeat k 0) a1)
                      (sub-ge (elimination-step A) (inc k)))))]
    (sub-ge A 0))))

(defn gaussian-elimination-with-pivoting 
  ([A b]
    (gaussian-elimination-with-pivoting (map #(conj (vec %1) %2) A b)))
  ([A]
    (letfn [(elimination-step [A]
              (map 
                (fn [ai]
                  (map 
                    #(- %1 (/ (* %2 (first ai)) (first (first A))))
                    (rest ai) (rest (first A)))) 
                (rest A)))
            (sub-ge [A k]
              (when-first [a1 A]
                (cons (concat (repeat k 0) a1)
                      (sub-ge (elimination-step ((pivot A) :A)) (inc k)))))]
    (sub-ge A 0))))

(defn lu-factorization-sparse [A]
  (letfn [(elimination-step [[A1 & Ak]]
            (reduce 
              (fn [AL [ai1 & aik]]
                (let [l (/ ai1 (first A1))
                      ai (mapv #(- %1 (* %2 l)) aik (rest A1))]
                  {:A (conj (AL :A) ai) :L (conj (AL :L) l)}))
              {:A [] :L []} Ak))
          (rearrange-L [L p]
            (let [pn (rest (reductions
                      (fn [[p1 & pk] i] (mapv #(if (< p1 %) (dec %) %) pk))
                      p p))]
              (mapv replace L pn)))]
    (loop [A A i (range (count A)) L [] U [] p []]
      (if (first A)
        (let [Ap (pivot A i 0)
              Al (elimination-step (Ap :A))
              lj (Al :L)
              ui (first (Ap :A))
              pi (first (Ap :p))]
          (recur (Al :A) (rest (Ap :p)) (conj L lj) (conj U ui) (conj p pi)))
      (zipmap [:L :U :p] [(rearrange-L L p) U p])))))


(defn lu-forward-substitution-sparse
  ([LUP b] (lu-forward-substitution-sparse (LUP :L) (LUP :p) b))
  ([L p b]
    (if (coll? (first b))
      ; If b is a matrix (transposed) calculate each column separately.
      (map (partial lu-forward-substitution-sparse L p) b)
      (let [b (map (partial nth b) p)]
        (reduce
          (fn [b k]
            (concat
              (take k b)
              (map 
                #(- %1 (* %2 (nth b (dec k)))) 
                (drop k b) (nth L (dec k)))))
          b (range 1 (count b)))))))

(defn inverse-matrix [A]
  (let [LUP (lu-factorization-sparse A)]   ; Calculate LU factorization
    (->>
      (diag (repeat (count A) 1))          ; Create unit matrix
      (lu-forward-substitution-sparse LUP) ; Evaluate forward substition
      (backward-substitution-sparse LUP)   ; Evaluate backward substition
      (apply map vector))))                ; Take transpose


(defn cholesky-decomposition [A]
  ; A must be symmetric posite-definite matrix.
  (reduce
    (fn [C ai]
      (let [lower-units
              (reduce 
                (fn [ci cj]
                  (conj ci (/ (- (nth ai (count ci))
                                 (apply + (map * ci cj)))
                              (last cj))))
                [] C)
            diag-unit
              (Math/sqrt (- (last ai) 
                            (apply + (map #(* % %) lower-units))))]
        (conj C (conj lower-units diag-unit))))
    [] (map take (range 1 (inc (count A))) A)))
