;(def A [[1 2 3 4][6 5 4 1][9 8 7 1][9 2 1 3]])
;(def A [[1 2 3][4 5 1][3 2 1]]) (def b [8 5 2])
(def A [[1.0E-32 1][1 1]]) (def b [1 2])

(defn pivot [A k]
  (reduce 
    #(if (< (second %2) (nth (first %1) k))
      %1 
      (cons 
        (nth %1 (first %2))
        (concat (take (first %2) %1) (drop (inc (first %2)) %1))))
    A (map vector (range) (map #(nth % k) A))))
(println (pivot A 1))

(defn pivot2 [A k]
  (let [pivot-row (key (apply max-key val (zipmap (range) (map #(nth % k) A))))]
    (concat 
      [(nth A pivot-row)]
      (take pivot-row A)
      (drop (inc pivot-row) A))))
(println (pivot2 A 1))

;(defn gaussian-elimination [A b]
;  (let [n (count A)]
;    (loop [A A b b k 0]
;      (if (= k (dec n)) 
;        [(triu A) b]
;        (let [a (map #(/ % (nth k (nth k A))) (map #(nth % k) (drop (inc k) A)))]
;          (recur (concat 
;                    (take k A) 
;                    (map #(drop (inc k) %) (drop (inc k) A)))
;                 b
;                 (inc k))))))

(defn gaussian-elimination [A b]
  (letfn [(A-ki-kj [A]
            (map 
              (fn [ai-jj]
                (map 
                  #(- %1 (/ (* %2 (first ai-jj)) (first (first A))))
                  (rest ai-jj) (rest (first A)))) 
              (rest A)))
          (recursive [A k]
            (when-first [a1 A]
              (cons (concat (repeat k 0) a1)
                    (recursive (A-ki-kj A) (inc k)))))]
  (recursive (map #(conj %1 %2) A b) 0)))
(println (gaussian-elimination A b))

(defn gaussian-elimination-with-pivoting [A b]
  (letfn [(A-ki-kj [A]
            (map 
              (fn [ai-jj]
                (map 
                  #(- %1 (/ (* %2 (first ai-jj)) (first (first A))))
                  (rest ai-jj) (rest (first A)))) 
              (rest A)))
          (recursive [A k]
            (when-first [a1 A]
              (cons (concat (repeat k 0) a1)
                    (recursive (A-ki-kj (pivot A 0)) (inc k)))))]
  (recursive (map #(conj %1 %2) A b) 0)))
(println (gaussian-elimination-with-pivoting A b))


